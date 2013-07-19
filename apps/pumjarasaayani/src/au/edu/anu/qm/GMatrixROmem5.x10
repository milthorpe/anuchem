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
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.util.ArrayList;
import x10.util.RailUtils;
import x10.util.Team;
import x10.util.WorkerLocalHandle;
import x10.lang.PlaceLocalHandle;
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
    val roK:Long; val roZ:Double; val omega:Double; val roThresh:Double;
    val place2atom:Rail[Long]; val place2func:Rail[Long]; val funcAtPlace:Rail[Long]; val offsetAtPlace:Rail[Long];

    // Standard conventional stuff
    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
    val nOrbital:Long; 
    val norm:Rail[Double]; 
    val emptyRailD = new Rail[Double](0), emptyRailI = new Rail[Long](0); 
    val shellPairs:PlaceLocalHandle[Rail[ShellPair]];
    val ylms:PlaceLocalHandle[Rail[Ylm]]; 
    val taux:WorkerLocalHandle[Integral_Pack];

    val place2ShellPair:Rail[Long];
    val numSigShellPairs:Long;
    val maxam1:Long;

    public def this(N:Long, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Long, omega:Double,roThresh:Double):GMatrixROmem5{self.M==N,self.N==N} {     
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
        } else { // full Coulomb operator
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
        place2atom = new Rail[Long](nPlaces+1);
        place2func = new Rail[Long](nPlaces+1);
        funcAtPlace = new Rail[Long](nPlaces);
        offsetAtPlace = new Rail[Long](nPlaces+1);
        var placeID:Long = nPlaces-1; var func:Long=0;
        val noOfAtoms=mol.getNumberOfAtoms();

        place2atom(nPlaces)=noOfAtoms-1;
        place2func(nPlaces)=mol.getAtom(noOfAtoms-1).getBasisFunctions().size();
        Console.OUT.printf("place %3d: atom=%5d function=%3d\n",nPlaces,place2atom(nPlaces),place2func(nPlaces));

        for(var a:Long=noOfAtoms-1; a>=0 && placeID>0; a--) { // centre a  
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();          
            for(var i:Long=naFunc-1; i>=0 && placeID>0; i--) { // basis functions on a
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
        var nShell:Long=0;
        for(var a:Long=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }
        Console.OUT.printf("nShell=%d...\n", nShell);

        val rawShellPairs = new Rail[ShellPair](nShell*nShell); // redundant list 
        val zeroPoint = Point3d(0.0, 0.0, 0.0);
       
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRailD, emptyRailD, emptyRailD, emptyRailD, 0, 0, 0, 0, emptyRailI, threshold);
        var mu:Long=0,nu:Long=0,ind:Long=0,totFunc:Long=0;
        val place2ShellPair=new Rail[Long](nPlaces+1);
        place2ShellPair(0)=0; placeID=1; 
        offsetAtPlace(0)=0; offsetAtPlace(nPlaces)=N;
        Console.OUT.printf("place %3d: offset = %5d shellpair #%5d\n",0,0,place2ShellPair(0));
        for(var a:Long=0; a<noOfAtoms; a++) { // centre a  
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();          
            for (var i:Long=0; i<naFunc; i++) { // basis functions on a
                val iaFunc = aFunc.get(i);    
                val aa=iaFunc.getTotalAngularMomentum();
                val aang = iaFunc.getTotalAngularMomentum();
                val aPoint = iaFunc.origin;
                val zetaA = iaFunc.exponents;
                val conA = iaFunc.coefficients; 
                val dConA = conA.size as Int;
                val maxbraa = (aa+1)*(aa+2)/2;
                if (a==place2atom(placeID) && i==place2func(placeID)) {
                    place2ShellPair(placeID)=ind;                    
                    funcAtPlace(placeID-1) = mu-offsetAtPlace(placeID-1);  
                    Console.OUT.printf("#basis functions1 = %5d #basis functions2 = %5d\n",funcAtPlace(placeID-1),totFunc); totFunc=0;
                    offsetAtPlace(placeID)= mu;          
                    Console.OUT.printf("place %3d: offset = %5d shellpair #%5d\n",placeID,offsetAtPlace(placeID),ind);
                    placeID++;
                }          
                for(var b:Long=0; b<noOfAtoms; b++) { // centre b
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = bFunc.size();                    
                    for(var j:Long=0; j<nbFunc; j++) { // basis functions on b
                        val jbFunc = bFunc.get(j);
                        val bb=jbFunc.getTotalAngularMomentum();                       
                        val maxbrab = (bb+1)*(bb+2)/2;     
                        val bang = jbFunc.getTotalAngularMomentum();
                        val bPoint = jbFunc.origin; 
                        val zetaB = jbFunc.exponents; 
                        val conB = jbFunc.coefficients; 
                        val dConB = conB.size as Int;
                        var contrib:Double = 0.; // ss = conservative estimate
                        val R2 = Math.pow(aPoint.i-bPoint.i,2.)+Math.pow(aPoint.j-bPoint.j,2.)+Math.pow(aPoint.k-bPoint.k,2.);
                        for (var ii:Int=0; ii<dConA; ii++) for (var jj:Int=0; jj<dConB; jj++) 
                            contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2)
                                     /Math.pow(jd.roZ,aang+bang);  // See Szabo Ostlund 3.284-3.286 
                        contrib=Math.abs(contrib);
                        if (contrib>=threshold) {
                            val maxL = new Rail[Long](roN+1); for (var ron:Long=0; ron<=roN; ron++) maxL(ron)=roL;
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
        Console.OUT.printf("#basis functions1 = %5d #basis functions2 = %5d\n", funcAtPlace(placeID-1), totFunc);
        Console.OUT.printf("place %3d: offset = %5d shellpair #%5d\n", nPlaces,mu+1, place2ShellPair(nPlaces));
      
        Console.OUT.printf("Found %d significant shellpairs.\n", numSigShellPairs);       
        
        val shellPairs = PlaceLocalHandle.make[Rail[ShellPair]](PlaceGroup.WORLD,
            ()=>new Rail[ShellPair](place2ShellPair(here.id+1)-place2ShellPair(here.id),
                (i:Long) => rawShellPairs(i+place2ShellPair(here.id))));
        this.place2ShellPair = place2ShellPair;
        this.shellPairs = shellPairs;

        Console.OUT.printf("Omega=%5.3f thresh=%e rad=%7.3f\n", omega, roThresh, jd.rad);
        val taux = new WorkerLocalHandle[Integral_Pack](() => new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        this.taux = taux;

        val roL_val = roL;
        val ylms = PlaceLocalHandle.make[Rail[Ylm]](PlaceGroup.WORLD, 
            ()=>new Rail[Ylm](place2ShellPair(here.id+1)-place2ShellPair(here.id),
                (i:Long) => 
                {
                    val sh = shellPairs()(i);
                    val tempY = new Rail[Double](sh.dconA*sh.dconB*(roL_val+1)*(roL_val+1));
                    taux().genClassY(sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.dconA, sh.dconB, roL_val, tempY);
                    new Ylm(tempY, roL_val)
                } ));
        this.ylms = ylms;

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
        val N=this.N;val nOrbital=this.nOrbital; val roN=this.roN; val roNK=this.roNK; val roL=this.roL; val roK=this.roK; val roZ=this.roZ; val norm=this.norm;
        val funcAtPlace=this.funcAtPlace; val offsetAtPlace=this.offsetAtPlace; val taux=this.taux;

        val TIMER_TOTAL = 0; val TIMER_JMATRIX = 1; val TIMER_KMATRIX = 2; val TIMER_GENCLASS = 3;
        
        timer.start(TIMER_TOTAL); 
        val jd = JobDefaults.getInstance();

        val maxTh=Runtime.NTHREADS; val nPlaces=Place.MAX_PLACES;            

        val cbs_auxInt= new Rail[Long](1);  cbs_auxInt(0)=N*roK;
        //@Ifdef("__DEBUG__") {  for (i in (0..(nPlaces-1))) Console.OUT.printf("%d %d\n",i,funcAtPlace(i)); Console.OUT.printf("%d\n",cbs_auxInt(0));}
        val auxIntGrid = new Grid(funcAtPlace, cbs_auxInt);
        val auxIntMat = DistDenseMatrix.make(auxIntGrid);

        val auxIntGrid2 = new Grid(cbs_auxInt,funcAtPlace);
        val auxIntMat2 = DistDenseMatrix.make(auxIntGrid2);

        val cbs_HalfAuxInt= new Rail[Long](1);  cbs_HalfAuxInt(0)=nOrbital*roK;
        val halfAuxGrid = new Grid(funcAtPlace, cbs_HalfAuxInt);
        val halfAuxMat = DistDenseMatrix.make(halfAuxGrid);

        val jMatrix = new DenseMatrix(N, N);
        val kMatrix = new DenseMatrix(N, N);
        val kval = GlobalRef(kMatrix);
        val jval = GlobalRef(jMatrix);
        val dk = new Rail[Double](roK); // eqn 15b in RO#7
        val dkval = GlobalRef(dk);
        val dummy = new Rail[Rail[Double]](0);
        val ttemp = new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](maxam1*maxam1*roK));
        val ttemp2= new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](maxam1*maxam1*roK));
        val tdk = new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](roK));
        val tjMatrix = new WorkerLocalHandle[DenseMatrix](() => new DenseMatrix(N,N));

        var tINT:Double=0.,tJ:Double=0.,tK:Double=0.;
        val gMat = new DenseMatrix(N,N);

        for (var ron:Long=0; ron<=roN; ron++)  {  
            // @Ifdef("__DEBUG__") { Console.OUT.printf("ron=%d\n",ron); }
                      
            timer.start(TIMER_GENCLASS); 
            dk.clear();
            val lron=ron; 

            // Distributed Generation of AuxMat
            // Console.OUT.println("Aux - distributed"); 
            finish ateach(place in Dist.makeUnique()) {
                val pid = here.id; val shp=shellPairs(); val ylmp = ylms();

                tdk.applyLocal((dk:Rail[Double]) => { dk.clear(); });
                val localMat=auxIntMat.local(); val localMat2=auxIntMat2.local();
                //@Ifdef("__DEBUG__") {Console.OUT.println("pid=" + pid + " starts..."); }
                finish for (thNo in 0..(maxTh-1)) async {
                    val aux = taux(); val temp=ttemp(); val temp2=ttemp2(); val myThreaddk=tdk();
                    
                    for (var spInd:Long=thNo+place2ShellPair(pid); spInd<place2ShellPair(pid+1); spInd+=maxTh) {
                        val sp=shp(spInd-place2ShellPair(pid)); val maxLron=sp.maxL(lron);                    
                        if (maxLron>=0) {

                            val maxLm=(maxLron+1)*(maxLron+1);
                            val y=ylmp(spInd-place2ShellPair(pid));
                            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, lron as Int, maxLron as Int, y.y, y.maxL); 


                            var ind:Long=0;
                            val musize=sp.mu2-sp.mu+1; val nusize=sp.nu2-sp.nu+1;
                            for (var tmu:Long=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) {
                                val scdmn=density(tmu,tnu); val nrm=norm(tmu)*norm(tnu); 
                                val ttmu=tmu-sp.mu; val ttnu=tnu-sp.nu; val tmuoff=tmu-offsetAtPlace(pid);
                                for (var rolm:Long=0; rolm<maxLm; rolm++) {
                                    val normAux = nrm*temp(ind++);       
                                    myThreaddk(rolm) += scdmn*normAux; 
                                    temp2((ttnu*roK+rolm)*musize+ttmu) = normAux;  
                                    localMat2(tnu*roK+rolm, tmuoff)=normAux;
                                }
                            }                    

                            ind=0;
                            for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) 
                                for (var rolm:Long=0; rolm<maxLm; rolm++) 
                                    for (var tmu:Long=sp.mu-offsetAtPlace(pid); tmu<=sp.mu2-offsetAtPlace(pid); tmu++)
                                        localMat(tmu,tnu*roK+rolm) = temp2(ind++);                            
                        }
                    }
                }

                //Console.OUT.println(pid + " Collecting dkvalue..."); 
                val srcDk = tdk.reduceLocal((a:Rail[Double],b:Rail[Double]) => RailUtils.map(a, b, a, (x:Double,y:Double)=>x+y));
                
                at(dkval.home) async {
                    val dstDk = (dkval());
                    atomic RailUtils.map(dstDk, srcDk, dstDk, (x:Double,y:Double)=>x+y);
                }
            }
            timer.stop(TIMER_GENCLASS); tINT+=(timer.last(TIMER_GENCLASS) as Double)/1e9;

            // J - distributed
            // Console.OUT.println("J - distributed"); 
            timer.start(TIMER_JMATRIX);   
            finish ateach(place in Dist.makeUnique()) {
                val pid = here.id; val shp=shellPairs();
                val localMat=auxIntMat2.local();
                tjMatrix.applyLocal((j:DenseMatrix) => { j.reset(); });
                finish for (thNo in 0..(maxTh-1)) async {
                    val myThreadJMat = tjMatrix();
                    for (var spInd:Long=thNo+place2ShellPair(pid); spInd<place2ShellPair(pid+1); spInd+=maxTh) {
                        val sp=shp(spInd-place2ShellPair(pid));
                        val maxLron=sp.maxL(lron);
                        if (sp.maxL(lron)>=0) { 
                            val maxLm=(maxLron+1)*(maxLron+1); 
                            for (var tmu:Long=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) {
                                var jContrib:Double=0.;  val tnuroK=tnu*roK; val tmuoff = tmu-offsetAtPlace(pid);
                                for (var rolm:Long=0; rolm<maxLm; rolm++) 
                                    jContrib += dk(rolm)*localMat(tnuroK+rolm, tmuoff);
                                myThreadJMat(tmu,tnu) += jContrib;
                            } 
                        }
                    }
                }

                val partialJMatrix = tjMatrix.reduceLocal((a:DenseMatrix,b:DenseMatrix) => a.cellAdd(b));

                at (jval.home) {
                    val jj=(jval());
                    // TODO is all of jMatrix required from each place?
                    for (var tnu:Long=0; tnu<N; tnu++) for (var tmu:Long=offsetAtPlace(pid); tmu<offsetAtPlace(pid+1); tmu++)
                    jj(tmu, tnu)+=partialJMatrix(tmu, tnu);
                }
                // Will convert to WorkerLocalHandle or PlaceLocalHandle
            }
            timer.stop(TIMER_JMATRIX); tJ+=(timer.last(TIMER_JMATRIX) as Double)/1e9;

            // K
            // Console.OUT.println("K - distributed"); 
            if (ron<=roNK) { // This produces K/2
                 timer.start(TIMER_KMATRIX);  



                 finish ateach(place in Dist.makeUnique())  {
                     val pid = here.id; //Console.OUT.println("pid=" + pid + " starts..."); 
                     val A=new DenseMatrix(funcAtPlace(pid)*roK, N, auxIntMat.local().d);
                     val B=new DenseMatrix(funcAtPlace(pid)*roK, nOrbital, halfAuxMat.local().d);
                     DenseMatrixBLAS.compMultTrans(A, mos, B, [funcAtPlace(pid)*roK, nOrbital, N], false);
                 }   

                 val mult=Math.ceil(nPlaces*.5+.5) as Long;
                 //Console.OUT.println("mult=" + mult);  
                 finish ateach(place in Dist.makeUnique())  {
                     val pid = here.id; //Console.OUT.println("pid=" + pid + " starts...");   
                     val a=halfAuxMat.local();                   
                     val moff=offsetAtPlace(pid);
                     for (var blk:Long=0; blk<mult; blk++) {
                         val b=at(Place((pid+blk)%nPlaces)) {halfAuxMat.local()};
                         val noff=offsetAtPlace((pid+blk)%nPlaces);
                         val c=new DenseMatrix(a.M,b.M);
                         c.multTrans(a, b, false);
                         //Console.OUT.printf("pid=%d, blk=%d [%d %d] x [%d %d]\n", pid, blk, a.M, a.N, b.M, b.N);
                         //Console.OUT.printf("moff=%d, noff=%d\n", moff, noff);
                         //c.debugPrint("c");
                         at (kval.home) {
                             for (var i:Long=0; i<c.M; i++) for (var j:Long=0; j<c.N; j++)
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
        /*finish ateach(place in Dist.makeUnique()) {
            val shp = shellPairs(); val numSp=place2ShellPair(here.id+1)-place2ShellPair(here.id);
            for (thNo in 0..(maxTh-1)) async for (var spInd:Long=thNo; spInd<numSp; spInd+=maxTh) {
                val sp=shp(spInd);
                if (sp.mu!=sp.nu) 
                at (jval.home) {
                    val jj=(jval());
                    for (var tmu:Long=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) 
                        jj(tnu,tmu) = jj(tmu,tnu);
                }
            }
        }*/

        val mult=Math.ceil(nPlaces*.5+.5);
        for (var p:Long=0; p<nPlaces; p++) for (var blk:Long=0; blk<mult; blk++) {
            val q=(p+blk)%nPlaces;
            for (var i:Long=offsetAtPlace(p); i<offsetAtPlace(p+1); i++) for (var j:Long=offsetAtPlace(q); j<offsetAtPlace(q+1); j++)                
               kMatrix(j,i) = kMatrix(i,j); 
        }


        // Use the upper half of J and K to form G
        finish for (thNo in 0..(maxTh-1)) async for (var tmu:Long=thNo; tmu<N; tmu+=maxTh) for (var tnu:Long=tmu; tnu<N; tnu++) 
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
            val eK = density.clone().mult(density, kMatrix).trace();
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