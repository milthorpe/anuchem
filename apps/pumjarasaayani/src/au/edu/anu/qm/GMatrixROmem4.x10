/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2013.
 */
package au.edu.anu.qm;

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;
import x10.compiler.Native;
import x10.compiler.NativeCPPInclude;
import x10.io.IOException;
import x10.util.ArrayList;
import x10.array.DistArray;
import x10.util.ArrayUtils;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;
import x10.matrix.dist.DistDenseMatrix;
import x10.matrix.dist.summa.SummaDense;

import x10x.vector.Vector;
import x10x.vector.Point3d;

import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ShellPair; 
import au.edu.anu.qm.Ylm; 
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;

import edu.utk.cs.papi.PAPI;

@NativeCPPInclude("mkl_math.h")
public class GMatrixROmem4 extends DenseMatrix{self.M==self.N} {

    // Timer & PAPI performance counters
    public val timer = new StatisticalTimer(4);

    // @Ifdef("__PAPI__") // XTENLANG-3132
    transient var papi:PAPI = new PAPI(); 

    // RO stuff 
    var roN:Int; var roNK:Int; var roL:Int; // 'var' because it can be overridden
    val roK:Int; val roZ:Double; val omega:Double; val roThresh:Double;
    val ylms:Rail[Ylm];   

    // Standard conventional stuff
    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
    val nOrbital:Int; val numSigShellPairs:Int;
    val norm:Rail[Double]; 
    val emptyRailD = new Rail[Double](0), emptyRailI = new Rail[Int](0); 
    val shellPairs:Rail[ShellPair]; 
    val place2ShellPair:Rail[Int];
    val maxam1:Int;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int, omega:Double,roThresh:Double):GMatrixROmem4{self.M==N,self.N==N} {     
        super(N, N);
        Console.OUT.printf("\nGMatrixROmem4.x10 'public def this' %s...\n", getDateString());
        this.bfs = bfs; this.mol = molecule; this.nOrbital = nOrbital; this.omega=omega; this.roThresh=roThresh;  
        // Set up RO N and L
        val jd = JobDefaults.getInstance();
        val l_n = new Rail[Int](jd.roN+3);
        val aux = new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ);
        if (omega>0.) { // long-range Ewald operator
            aux.getNL(l_n);
            roN=roNK=l_n(0);
            roL=l_n(roN+2);  
            this.roNK=roN;
        }
        else { // full Coulomb operator
            this.roN=jd.roN;
            this.roL=jd.roL;
            if (jd.roNK==-1) this.roNK=roN; else this.roNK=jd.roNK; 
        }
        roK = (roL+1)*(roL+1); roZ=jd.roZ;      

        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val mdc=bfs.getShellList().getMaximumDegreeOfContraction();
        maxam1 = (maxam+1)*(maxam+2)/2;
        this.norm = bfs.getNormalizationFactors(); // Vector

        // Shell/Shellpair business
        Console.OUT.printf("GMatrixROmem4.x10 Screening shellpairs...\n");
        val threshold = roThresh*jd.roZ*jd.roZ*1e-3; // ** must be relative to roThresh *** otherwise Z scaling will cause a problem
        // This is effectively a density threshold RO Thesis (2.26)
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }
        Console.OUT.printf("nShell=%d...\n", nShell);

        var rawShellPairs:Rail[ShellPair] = new Array[ShellPair](nShell*nShell); // redundant list 
        val zeroPoint = Point3d(0.0, 0.0, 0.0);
       
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRailD, emptyRailD, emptyRailD, emptyRailD, 0, 0, 0, 0, emptyRailI, threshold);
        var mu:Int=0,nu:Int=0,ind:Int=0,totFunc:Int=0;
        for(var a:Int=0; a<noOfAtoms; a++) { // centre a  
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();          
            for(var i:Int=0; i<naFunc; i++) { // basis functions on a
                val iaFunc = aFunc.get(i);               
                for(var b:Int=0; b<noOfAtoms; b++) { // centre b
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = bFunc.size();                    
                    for(var j:Int=0; j<nbFunc; j++) { // basis functions on b
                        val jbFunc = bFunc.get(j);
                        var aaFunc:ContractedGaussian=iaFunc, bbFunc:ContractedGaussian=jbFunc;
                        val aa=iaFunc.getTotalAngularMomentum(); val bb=jbFunc.getTotalAngularMomentum();                       
                        val maxbraa = (aa+1)*(aa+2)/2; val maxbrab = (bb+1)*(bb+2)/2;     
                        // Symetric with respect to mu/nu - but we make the list redundant 
                        // Asymmetric
                        if (mu>=nu) {
                        val aang = aaFunc.getTotalAngularMomentum(); val bang = bbFunc.getTotalAngularMomentum();
                        val aPoint = aaFunc.origin; val bPoint = bbFunc.origin; 
                        val zetaA = aaFunc.exponents; val zetaB = bbFunc.exponents; 
                        val conA = aaFunc.coefficients; val conB = bbFunc.coefficients; 
                        val dConA = conA.size; val dConB = conB.size;
                        var contrib : Double = 0.; // ss = conservative estimate
                        val R2 = Math.pow(aPoint.i-bPoint.i,2.)+Math.pow(aPoint.j-bPoint.j,2.)+Math.pow(aPoint.k-bPoint.k,2.);
                        for (var ii:Int=0; ii<dConA; ii++) for (var jj:Int=0; jj<dConB; jj++) 
                            contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2)
                                     /Math.pow(jd.roZ,aang+bang);  // See Szabo Ostlund 3.284-3.286 
                        contrib=Math.abs(contrib);
                        if (contrib>=threshold) {
                            val maxL = new Rail[Int](roN+1); for (var ron:Int=0; ron<=roN; ron++) maxL(ron)=roL;
                            rawShellPairs(ind) = new ShellPair(aang,bang,aPoint,bPoint,zetaA,zetaB,conA,conB,dConA,dConB,mu,nu,maxL,contrib);     
                            ind++;
                            totFunc+=maxbraa*maxbrab; // *(roN+1)*(roL+1)*(roL+1)
                        }               
                        } // Asym
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab; else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        val nPlaces:Int=Place.MAX_PLACES;
        val fpp:Int=Math.ceil(totFunc/nPlaces) as Int;
        place2ShellPair=new Rail[Int](nPlaces+1);
        var placeID:Int=nPlaces-1, func:Int=0;
        Console.OUT.printf("totFunc=%d, nPlace=%d, fpp=%d, ind=%d...\n", totFunc, nPlaces, fpp, ind);
        for (var spID:Int=ind-1; spID>0 && placeID>0; spID--) {
             val sp=rawShellPairs(spID);
             val aa=sp.aang, bb=sp.bang;
             val maxbraa = (aa+1)*(aa+2)/2; val maxbrab = (bb+1)*(bb+2)/2;
             func+=maxbraa*maxbrab;
             if (totFunc-func<placeID*fpp) {
                 Console.OUT.printf("place %d : spID %d...\n", placeID, spID);
                 place2ShellPair(placeID--)=spID;
             }
        }
        place2ShellPair(0)=0; place2ShellPair(nPlaces)=ind;
        Console.OUT.printf("***By default place 0 : spID 0***\n");
        // if there are too few shellpairs this might break down
        // Should check integrity of the list  place2ShellPair

        this.numSigShellPairs=ind;
        Console.OUT.printf("Found %d significant shellpairs.\n",numSigShellPairs);
        this.shellPairs = new Rail[ShellPair](numSigShellPairs);    
        // this calculation should be distributed too
        this.ylms = new Rail[Ylm](numSigShellPairs);     
        for (i in 0..(numSigShellPairs-1)) {
            val sh=rawShellPairs(i); 
            this.shellPairs(i)=sh;            
            val tempY=new Rail[Double](sh.dconA*sh.dconB*(roL+1)*(roL+1));
            aux.genClassY(sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.dconA, sh.dconB, roL, tempY);
            ylms(i) = new Ylm(tempY,roL);
        }
        rawShellPairs = null; // Deallocate this variable

        @Ifdef("__MKL__") {
            Console.OUT.print("mklGetMaxThreads() was " + mklGetMaxThreads() + " and is now set to"); mklSetNumThreads(Runtime.NTHREADS);
            Console.OUT.println(" " + mklGetMaxThreads() + " thread(s).");
        }  

        @Ifdef("__PAPI__") { 
            papi.initialize();
            papi.countFlops();
            papi.countMemoryOps(); 
        }
    }

    @Native("c++", "mkl_get_max_threads()") private native static def mklGetMaxThreads():Int;
    @Native("c++", "MKL_Set_Num_Threads(#a)") private native static def mklSetNumThreads(a:Int):void;

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        Console.OUT.printf("\nGMatrixROmem.x10 'public def compute' %s...\n", getDateString()); 
        val timer=this.timer; val shellPairs=this.shellPairs; val maxam1=this.maxam1; val numSigShellPairs=this.numSigShellPairs; val ylms=this.ylms;
        val N=this.N;val nOrbital=this.nOrbital; val roN=this.roN; val roNK=this.roNK; val roL=this.roL; val roK=this.roK; val roZ=this.roZ; val omega=this.omega; val roThresh=this.roThresh; val norm=this.norm;
        val place2ShellPair=this.place2ShellPair;
        val TIMER_TOTAL = 0; val TIMER_JMATRIX = 1; val TIMER_KMATRIX = 2; val TIMER_GENCLASS = 3;
        
        timer.start(TIMER_TOTAL); 
        val jd = JobDefaults.getInstance();
        this.reset(); val gVal = GlobalRef(this); 

        val maxTh=Runtime.NTHREADS; val maxPl=Place.MAX_PLACES;            

        val auxIntMat = DistDenseMatrix.make(N*roK,N);
        val halfAuxMat = DistDenseMatrix.make(nOrbital,N*roK);
        val jMatrix = new DenseMatrix(N, N);
        val kMatrix = DistDenseMatrix.make(N, N);
        val dk = new Rail[Double](roK); // eqn 15b in RO#7
        val dkval =GlobalRef(dk);
        val ttemp = new Rail[Rail[Double]](maxTh, (Int) => new Rail[Double](maxam1*maxam1*roK));
        val tdk = new Rail[Rail[Double]](maxTh, (Int) => new Rail[Double](roK));
        val tjMatrix = new Rail[DenseMatrix](maxTh, (Int) => new DenseMatrix(N,N));

        var tINT:Double=0.,tJ:Double=0.,tK:Double=0.;
        val gMat = new DenseMatrix(N,N);

        for (var ron:Int=0; ron<=roN; ron++)  {  
            @Ifdef("__DEBUG__") {Console.OUT.printf("ron=%d...\n",ron); }
            finish for (thNo in 0..(maxTh-1)) async tdk(thNo).clear();     
                      
            timer.start(TIMER_GENCLASS); 
dk.clear();
val lron=ron; 
// Local variable is accessed at a different place, and therefore it must be an initialized val.       [exec]      	 Variable name: ron
            finish ateach(place in Dist.makeUnique()) async {
                val pid = here.id; Console.OUT.println("pid=" + pid); 

                val taux = new Rail[Integral_Pack](maxTh, (Int) => new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
                finish for (thNo in 0..(maxTh-1)) async {
                    val aux = taux(thNo); 
                    for (var spInd:Int=thNo+place2ShellPair(pid); spInd<place2ShellPair(pid+1); spInd+=maxTh) {
                        val sp=shellPairs(spInd); val maxLron=sp.maxL(lron);                    
                        if (maxLron>=0) {
                            val maxLm=(maxLron+1)*(maxLron+1); var ind:Int=0; val temp=ttemp(thNo); val myThreaddk=tdk(thNo);  
                            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, lron, maxLron,ylms(spInd).y, ylms(spInd).maxL); 
                            if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) {
                                val scdmn=density(tmu,tnu)*2.; val nrm=norm(tmu)*norm(tnu); 
                                val tmuroK=tmu*roK; val tnuroK=tnu*roK;
                                for (var rolm:Int=0; rolm<maxLm; rolm++) {
                                    val normAux = nrm*temp(ind++); 
                                    myThreaddk(rolm) += scdmn*normAux; 
                                    auxIntMat(tmuroK+rolm,tnu) = auxIntMat(tnuroK+rolm,tmu) = normAux;                            
                                } 
                            }
                            else for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) {
                                val scdmn=density(tmu,tnu); val nrm=norm(tmu)*norm(tnu);
                                val tmuroK=tmu*roK; 
                                for (var rolm:Int=0; rolm<maxLm; rolm++) {
                                    val normAux = nrm*temp(ind++);
                                    myThreaddk(rolm) += scdmn*normAux; 
                                    auxIntMat(tmuroK+rolm,tnu) = normAux;                            
                                } 
                            }
                        }
                    }
                }
                
                finish for (thNo in 0..(maxTh-1)) async for (thNo2 in 0..(maxTh-1)) {
                    val myThreaddk = tdk(thNo2);
                    for (var rolm:Int=thNo; rolm<roK; rolm+=maxTh) { val rolmm=rolm;
                        //at (gVal) async atomic dk(rolmm)+=myThreaddk(rolmm);}
                        at(dkval.home) async atomic (dkval())(rolmm)+=myThreaddk(rolmm); }
                }
            }
            timer.stop(TIMER_GENCLASS); tINT+=(timer.last(TIMER_GENCLASS) as Double)/1e9;

            // J
            timer.start(TIMER_JMATRIX);   
            // add for at async    
            finish for (thNo in 0..(maxTh-1)) async tjMatrix(thNo).reset();
            finish for (thNo in 0..(maxTh-1)) async {
                val myThreadJMat = tjMatrix(thNo);      
                for (var spInd:Int=thNo; spInd<numSigShellPairs; spInd+=maxTh) {
                    val sp=shellPairs(spInd);
                    val maxLron=sp.maxL(ron);
                    if (sp.maxL(ron)>=0) { 
                        val maxLm=(maxLron+1)*(maxLron+1); 
                        for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) {
                            var jContrib:Double=0.; val tmuroK=tmu*roK;
                            for (var rolm:Int=0; rolm<maxLm; rolm++) jContrib+=dk(rolm)*auxIntMat(tmuroK+rolm,tnu);
                                myThreadJMat(tmu,tnu) += jContrib;
                        } 
                    }
                }
            }

            finish for (thNo in 0..(maxTh-1)) async for (thNo2 in 0..(maxTh-1)) {
                val myThreadJMat = tjMatrix(thNo2);
                for (var tnu:Int=thNo; tnu<N; tnu+=maxTh) for (var tmu:Int=0; tmu<N; tmu++)
                    jMatrix(tmu,tnu)+=myThreadJMat(tmu,tnu);
            }

            timer.stop(TIMER_JMATRIX); tJ+=(timer.last(TIMER_JMATRIX) as Double)/1e9;

            // K - change to SUMMA
            if (ron<=roNK) { //This produces K/2
                 timer.start(TIMER_KMATRIX);  
                 //DenseMatrixBLAS.compMultTrans(mos, auxIntMat, halfAuxMat, [nOrbital,N*roK, N], false);

                 //This may present futhur difficulty when we change to DistMatrix
                 //val halfAuxMat2 = new DenseMatrix(roK*nOrbital, N, halfAuxMat.d);

                 //DenseMatrixBLAS.compTransMult(halfAuxMat2, halfAuxMat2, kMatrix, [N, N, roK*nOrbital], true);
                 timer.stop(TIMER_KMATRIX); tK+=(timer.last(TIMER_KMATRIX) as Double)/1e9;
            }
        
            // Fix upper half of J (see the definition of shellPairs) ==> sp.mu>sp.nu
            // Fix the whole matrix ==> if (sp.mu!=sp.nu)
            finish for (thNo in 0..(maxTh-1)) async for (var spInd:Int=thNo; spInd<numSigShellPairs; spInd+=maxTh) {
                val sp=shellPairs(spInd);
                if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) 
                    jMatrix(tnu,tmu) = jMatrix(tmu,tnu);
            }
            // Use the upper half of J and K to form G
            finish for (thNo in 0..(maxTh-1)) async for (var tmu:Int=thNo; tmu<N; tmu+=maxTh) for (var tnu:Int=tmu; tnu<N; tnu++) 
                gMat(tnu,tmu)+=gMat(tmu,tnu)=jMatrix(tmu,tnu)-kMatrix(tmu,tnu);
            
            //at(gVal) async atomic gVal().cellAdd(gMat); 
            //at(gVal) async atomic gVal().cellAdd(jMatrix); 
            Console.OUT.printf("Time INT = %.2f s J = %.2f s K = %.2f s\n", tINT, tJ, tK);
            Console.OUT.flush();
        }

        @Ifdef("__DEBUG__") {
        val eJ = density.clone().mult(density, jMatrix).trace();
        Console.OUT.printf("  EJ = %.10f a.u.\n", .25*eJ/jd.roZ);
        }

        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }

        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);
    }

    private def getDateString() {
        val result:String;
        try {
            val dateReader = Runtime.execForRead("date"); 
            result = dateReader.readLine();
        } catch (e:IOException) {
            // could not read date! use current time in milliseconds
            result = System.currentTimeMillis() + "ms";
        }
        return result;
    }
}
