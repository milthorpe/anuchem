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
import x10.matrix.block.Grid;
import x10.matrix.block.MatrixBlock;
import x10.matrix.block.DenseBlock;
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
public class GMatrixROmem5 extends DenseMatrix{self.M==self.N} {

    // Timer & PAPI performance counters
    public val timer = new StatisticalTimer(4);

    // @Ifdef("__PAPI__") // XTENLANG-3132
    transient var papi:PAPI = new PAPI(); 

    // RO stuff 
    var roN:Int; var roNK:Int; var roL:Int; // 'var' because it can be overridden
    val roK:Int; val roZ:Double; val omega:Double; val roThresh:Double;
    val ylms:Rail[Ylm];   
    val place2atom:Rail[Int]; val place2func:Rail[Int]; val funcAtPlace:Rail[Int]; val offsetAtPlace:Rail[Int];

    // Standard conventional stuff
    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
    val nOrbital:Int; 
    val norm:Rail[Double]; 
    val emptyRailD = new Rail[Double](0), emptyRailI = new Rail[Int](0); 
    val shellPairs:Rail[ShellPair]; 
    val place2ShellPair:Rail[Int];
    val numSigShellPairs:Int;
    val maxam1:Int;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int, omega:Double,roThresh:Double):GMatrixROmem5{self.M==N,self.N==N} {     
        super(N, N);
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def this' %s...\n", getDateString());
        this.bfs = bfs; this.mol = molecule; this.nOrbital = nOrbital; this.omega=omega; this.roThresh=roThresh;  

        // Set up RO N, L, K, Z
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

        // Some other handy variables
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val mdc=bfs.getShellList().getMaximumDegreeOfContraction();
        maxam1 = (maxam+1)*(maxam+2)/2;
        this.norm = bfs.getNormalizationFactors(); // Vector

        // Data distribution on a shell basis
        // Input: nPlaces, mol
        // Output: place2atom, place2func
        val nPlaces=Place.MAX_PLACES; 
        val npp=N/nPlaces;
        place2atom = new Rail[Int](nPlaces+1);
        place2func = new Rail[Int](nPlaces+1);
        funcAtPlace = new Rail[Int](nPlaces);
        offsetAtPlace = new Rail[Int](nPlaces+1);
        var placeID:Int = nPlaces-1; var func:Int=0;
        val noOfAtoms=mol.getNumberOfAtoms();

        place2atom(nPlaces)=noOfAtoms-1;
        place2func(nPlaces)=mol.getAtom(noOfAtoms-1).getBasisFunctions().size();
        Console.OUT.printf("place %3d: atom=%5d function=%3d\n",nPlaces,place2atom(nPlaces),place2func(nPlaces));

        for(var a:Int=noOfAtoms-1; a>=0 && placeID>0; a--) { // centre a  
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();          
            for(var i:Int=naFunc-1; i>=0 && placeID>0; i--) { // basis functions on a
                val iaFunc = aFunc.get(i);   
                val aa=iaFunc.getTotalAngularMomentum();
                func+=(aa+1)*(aa+2)/2;
                if (N-func<placeID*npp) {
                     place2atom(placeID)=a;
                     place2func(placeID)=i;
                     Console.OUT.printf("place %3d: atom=%5d function=%3d\n",placeID,a,i);
                     placeID--;
                }
            }
        }

        place2atom(0)=0;
        place2func(0)=0;
        Console.OUT.printf("place %3d: atom=%5d function=%3d\n",0,place2atom(0),place2func(0));

        // Generate shellpair
        val threshold = roThresh*jd.roZ*jd.roZ*1e-3; 
        // ** must be relative to roThresh *** otherwise Z scaling will cause a problem
        // This is effectively a density threshold RO Thesis (2.26)
        var nShell:Int=0;
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }
        Console.OUT.printf("nShell=%d...\n", nShell);

        var rawShellPairs:Rail[ShellPair] = new Array[ShellPair](nShell*nShell); // redundant list 
        val zeroPoint = Point3d(0.0, 0.0, 0.0);
       
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRailD, emptyRailD, emptyRailD, emptyRailD, 0, 0, 0, 0, emptyRailI, threshold);
        var mu:Int=0,nu:Int=0,ind:Int=0,totFunc:Int=0;
        place2ShellPair=new Rail[Int](nPlaces+1);
        place2ShellPair(0)=0; placeID=1; 
        offsetAtPlace(0)=0; offsetAtPlace(nPlaces)=N;
        Console.OUT.printf("place %3d: offset=%5d shellpair #%5d\n",0,0,place2ShellPair(0));
        for(var a:Int=0; a<noOfAtoms; a++) { // centre a  
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();          
            for (var i:Int=0; i<naFunc; i++) { // basis functions on a
                val iaFunc = aFunc.get(i);    
                val aa=iaFunc.getTotalAngularMomentum();
                val aang = iaFunc.getTotalAngularMomentum();
                val aPoint = iaFunc.origin;
                val zetaA = iaFunc.exponents;
                val conA = iaFunc.coefficients; 
                val dConA = conA.size;
                val maxbraa = (aa+1)*(aa+2)/2;
                if (a==place2atom(placeID) && i==place2func(placeID)) {
                    place2ShellPair(placeID)=ind;                    
                    funcAtPlace(placeID-1) = mu-offsetAtPlace(placeID-1);  
                    Console.OUT.printf("#basis fucntions1 =%5d #basis fucntions2 =%5d\n",funcAtPlace(placeID-1),totFunc); totFunc=0;
                    offsetAtPlace(placeID)= mu;          
                    Console.OUT.printf("place %3d: offset=%5d shellpair #%5d\n",placeID,offsetAtPlace(placeID),ind);
                    placeID++;
                }          
                for(var b:Int=0; b<noOfAtoms; b++) { // centre b
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = bFunc.size();                    
                    for(var j:Int=0; j<nbFunc; j++) { // basis functions on b
                        val jbFunc = bFunc.get(j);
                        val bb=jbFunc.getTotalAngularMomentum();                       
                        val maxbrab = (bb+1)*(bb+2)/2;     
                        val bang = jbFunc.getTotalAngularMomentum();
                        val bPoint = jbFunc.origin; 
                        val zetaB = jbFunc.exponents; 
                        val conB = jbFunc.coefficients; 
                        val dConB = conB.size;
                        var contrib:Double = 0.; // ss = conservative estimate
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
                        //} // Asym
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab; else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        numSigShellPairs=place2ShellPair(placeID)=ind; 
        funcAtPlace(placeID-1) = mu-offsetAtPlace(placeID-1); 
        Console.OUT.printf("#basis fucntions1 =%5d #basis fucntions2 =%5d\n", funcAtPlace(placeID-1), totFunc);
        Console.OUT.printf("place %3d: offset=%5d shellpair #%5d\n", nPlaces,mu+1, place2ShellPair(nPlaces));
      
        Console.OUT.printf("Found %d significant shellpairs.\n", numSigShellPairs);
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
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def compute' %s...\n", getDateString()); 
        val timer=this.timer; 
        val place2ShellPair=this.place2ShellPair; val shellPairs=this.shellPairs; val maxam1=this.maxam1; val numSigShellPairs=this.numSigShellPairs; val ylms=this.ylms;
        val N=this.N;val nOrbital=this.nOrbital; val roN=this.roN; val roNK=this.roNK; val roL=this.roL; val roK=this.roK; val roZ=this.roZ; val omega=this.omega; val roThresh=this.roThresh; val norm=this.norm;
        val funcAtPlace=this.funcAtPlace; val offsetAtPlace=this.offsetAtPlace;
        val TIMER_TOTAL = 0; val TIMER_JMATRIX = 1; val TIMER_KMATRIX = 2; val TIMER_GENCLASS = 3;
        
        timer.start(TIMER_TOTAL); 
        val jd = JobDefaults.getInstance();
        //this.reset(); val gVal = GlobalRef(this); 

        val maxTh=Runtime.NTHREADS; val nPlaces=Place.MAX_PLACES;            

        val cbs_auxInt= new Rail[Int](1);  cbs_auxInt(0)=N*roK;
        //@Ifdef("__DEBUG__") {  for (i in (0..(nPlaces-1))) Console.OUT.printf("%d %d\n",i,funcAtPlace(i)); Console.OUT.printf("%d\n",cbs_auxInt(0));}
        val auxIntGrid = new Grid(funcAtPlace, cbs_auxInt);
        val auxIntMat = DistDenseMatrix.make(auxIntGrid);

        val jMatrix = new DenseMatrix(N, N);
        val kMatrix = new DenseMatrix(N, N);
        val kval = GlobalRef(kMatrix);
        val jval = GlobalRef(jMatrix);
        val dk = new Rail[Double](roK); // eqn 15b in RO#7
        val dkval = GlobalRef(dk);
        val ttemp = new Rail[Rail[Double]](maxTh, (Int) => new Rail[Double](maxam1*maxam1*roK));
        val tdk = new Rail[Rail[Double]](maxTh, (Int) => new Rail[Double](roK));
        val tjMatrix = new Rail[DenseMatrix](maxTh, (Int) => new DenseMatrix(N,N));

        var tINT:Double=0.,tJ:Double=0.,tK:Double=0.;
        val gMat = new DenseMatrix(N,N);

        for (var ron:Int=0; ron<=roN; ron++)  {  
            // @Ifdef("__DEBUG__") { Console.OUT.printf("ron=%d\n",ron); }
            finish for (thNo in 0..(maxTh-1)) async tdk(thNo).clear();     
                      
            timer.start(TIMER_GENCLASS); 
            dk.clear();
            val lron=ron; 

            // Distributed Generation of AuxMat
            // Console.OUT.println("Aux - distributed"); 
            finish ateach(place in Dist.makeUnique()) {
                val pid = here.id; 
                val localMat=auxIntMat.local();
                //@Ifdef("__DEBUG__") {Console.OUT.println("pid=" + pid + " starts..."); }
                val taux = new Rail[Integral_Pack](maxTh, (Int) => new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
                finish for (thNo in 0..(maxTh-1)) async {
                    val aux = taux(thNo); 
                    for (var spInd:Int=thNo+place2ShellPair(pid); spInd<place2ShellPair(pid+1); spInd+=maxTh) {
                        val sp=shellPairs(spInd); val maxLron=sp.maxL(lron);                    
                        if (maxLron>=0) {
                            val maxLm=(maxLron+1)*(maxLron+1); var ind:Int=0; val temp=ttemp(thNo); val myThreaddk=tdk(thNo);  
                            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, lron, maxLron,ylms(spInd).y, ylms(spInd).maxL); 
                            for (var tmu:Int=sp.mu-offsetAtPlace(pid); tmu<=sp.mu2-offsetAtPlace(pid); tmu++) 
                            for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) {
                                val scdmn=density(tmu,tnu) ; val nrm=norm(tmu)*norm(tnu); 
                                for (var rolm:Int=0; rolm<maxLm; rolm++) {
                                    val normAux = nrm*temp(ind++); 
                                    myThreaddk(rolm) += scdmn*normAux; 
                                    localMat(tmu,tnu*roK+rolm) = normAux;
                                } 
                            }
                        }
                    }
                }

                //Console.OUT.println(pid + " Collecting dkvalue..."); 
                finish for (thNo in 0..(maxTh-1)) async for (thNo2 in 0..(maxTh-1)) {
                    val myThreaddk = tdk(thNo2);
                    for (var rolm:Int=thNo; rolm<roK; rolm+=maxTh) { val rolmm=rolm;
                        //at (gVal) async atomic dk(rolmm)+=myThreaddk(rolmm);}
                        at(dkval.home) async atomic (dkval())(rolmm)+=myThreaddk(rolmm); 
                    }
                }
            }
            timer.stop(TIMER_GENCLASS); tINT+=(timer.last(TIMER_GENCLASS) as Double)/1e9;

            // J - distributed
            // Console.OUT.println("J - distributed"); 
            timer.start(TIMER_JMATRIX);   
            finish ateach(place in Dist.makeUnique()) {
                val pid = here.id; 
                val localMat=auxIntMat.local();      
                finish for (thNo in 0..(maxTh-1)) async tjMatrix(thNo).reset();
                finish for (thNo in 0..(maxTh-1)) async {
                    val myThreadJMat = tjMatrix(thNo);      
                    for (var spInd:Int=thNo+place2ShellPair(pid); spInd<place2ShellPair(pid+1); spInd+=maxTh) {
                        val sp=shellPairs(spInd);
                        val maxLron=sp.maxL(lron);
                        if (sp.maxL(lron)>=0) { 
                            val maxLm=(maxLron+1)*(maxLron+1); 
                            for (var tmu:Int=sp.mu-offsetAtPlace(pid); tmu<=sp.mu2-offsetAtPlace(pid); tmu++) 
                            for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) {
                                var jContrib:Double=0.;  val tnuroK=tnu*roK;
                                for (var rolm:Int=0; rolm<maxLm; rolm++) 
                                    jContrib += dk(rolm)*localMat(tmu, tnuroK+rolm);
                                myThreadJMat(tmu,tnu) += jContrib;
                            } 
                        }
                    }
                }

                for (thNo in 1..(maxTh-1)) {
                    val myThreadJMat = tjMatrix(thNo);
                    for (var tnu:Int=0; tnu<N; tnu++) for (var tmu:Int=offsetAtPlace(pid); tmu<offsetAtPlace(pid+1); tmu++)
                        (tjMatrix(0))(tmu,tnu)+=myThreadJMat(tmu, tnu);
                }
                val myThreadJMat = tjMatrix(0);
                at (jval.home) for (var tnu:Int=0; tnu<N; tnu++) for (var tmu:Int=offsetAtPlace(pid); tmu<offsetAtPlace(pid+1); tmu++)
                    (jval())(tmu, tnu)+=myThreadJMat(tmu, tnu);

            }
            timer.stop(TIMER_JMATRIX); tJ+=(timer.last(TIMER_JMATRIX) as Double)/1e9;

            // K
            // Console.OUT.println("K - distributed"); 
            if (ron<=roNK) { // This produces K/2
                 timer.start(TIMER_KMATRIX);  

                 val cbs_HalfAuxInt= new Rail[Int](1);  cbs_HalfAuxInt(0)=nOrbital*roK;
                 val halfAuxGrid = new Grid(funcAtPlace, cbs_HalfAuxInt);
                 val halfAuxMat = DistDenseMatrix.make(halfAuxGrid);

                 finish ateach(place in Dist.makeUnique())  {
                     val pid = here.id; //Console.OUT.println("pid=" + pid + " starts..."); 
                     val A=new DenseMatrix(funcAtPlace(pid)*roK, N, auxIntMat.local().d);
                     val B=new DenseMatrix(funcAtPlace(pid)*roK, nOrbital, halfAuxMat.local().d);
                     DenseMatrixBLAS.compMultTrans(A, mos, B, [funcAtPlace(pid)*roK, nOrbital, N], false);
                 }   

                 val mult=Math.ceil(nPlaces*.5+.5) as Int;
                 //Console.OUT.println("mult=" + mult);  
                 finish ateach(place in Dist.makeUnique())  {
                     val pid = here.id; //Console.OUT.println("pid=" + pid + " starts...");   
                     val a=halfAuxMat.local();                   
                     val moff=offsetAtPlace(pid);
                     for (var blk:Int=0; blk<mult; blk++) {
                         val b=at(Place((pid+blk)%nPlaces)) {halfAuxMat.local()};
                         val noff=offsetAtPlace((pid+blk)%nPlaces);
                         val c=new DenseMatrix(a.M,b.M);
                         c.multTrans(a, b, false);
                         //Console.OUT.printf("pid=%d, blk=%d [%d %d] x [%d %d]\n", pid, blk, a.M, a.N, b.M, b.N);
                         //Console.OUT.printf("moff=%d, noff=%d\n", moff, noff);
                         //c.debugPrint("c");
                         at (kval.home) {
                             for (var i:Int=0; i<c.M; i++) for (var j:Int=0; j<c.N; j++)
                                 (kval())(moff+i,noff+j)+=c(i,j);
                         }
                     }
                 }
          
                 timer.stop(TIMER_KMATRIX); tK+=(timer.last(TIMER_KMATRIX) as Double)/1e9;
            }
        }

        // G - at place 0 - muti-threading
        Console.OUT.printf("\nG matrix\n"); 
        // Fix upper half of J (see the definition of shellPairs) ==> sp.mu>sp.nu
        // Fix the whole matrix ==> if (sp.mu!=sp.nu)
        finish for (thNo in 0..(maxTh-1)) async for (var spInd:Int=thNo; spInd<numSigShellPairs; spInd+=maxTh) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) 
                jMatrix(tnu,tmu) = jMatrix(tmu,tnu);
        }

        val mult=Math.ceil(nPlaces*.5+.5);
        for (var p:Int=0; p<nPlaces; p++) for (var blk:Int=0; blk<mult; blk++) {
            val q=(p+blk)%nPlaces;
            for (var i:Int=offsetAtPlace(p); i<offsetAtPlace(p+1); i++) for (var j:Int=offsetAtPlace(q); j<offsetAtPlace(q+1); j++)                
               kMatrix(j,i) = kMatrix(i,j); 
        }


        // Use the upper half of J and K to form G
        finish for (thNo in 0..(maxTh-1)) async for (var tmu:Int=thNo; tmu<N; tmu+=maxTh) for (var tnu:Int=tmu; tnu<N; tnu++) 
            gMat(tnu,tmu)=gMat(tmu,tnu)=jMatrix(tmu,tnu)-kMatrix(tmu,tnu);

        this.reset();
        this.cellAdd(gMat);

        Console.OUT.printf("Time INT = %.2f s J = %.2f s K = %.2f s\n", tINT, tJ, tK);
        Console.OUT.flush();

        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);
 
        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }

        @Ifdef("__DEBUG__") {
            val eJ = density.clone().mult(density, jMatrix).trace();
            val kDmat = new DenseMatrix(N,N);
            for ([i,j] in ( (0..(N-1))*(0..(N-1)) ) )
                kDmat(i,j)=kMatrix(i,j);
            val eK = density.clone().mult(density, kDmat).trace();
            Console.OUT.printf("  EJ = %.10f EK=%.10f\n", .25*eJ/jd.roZ, .25*eK/jd.roZ);
        }       
 
        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }

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
