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
import x10.compiler.Native;
import x10.compiler.NativeCPPInclude;
import x10.io.IOException;

import x10.util.RailUtils;
import x10.util.Team;

import x10.util.WorkerLocalHandle;
import x10.lang.PlaceLocalHandle;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;
import x10.matrix.dist.DistDenseMatrix;
import x10.matrix.block.Grid;
import x10.regionarray.Dist;

import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ShellPair; 
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;

import edu.utk.cs.papi.PAPI;

@NativeCPPInclude("mkl_math.h")

public class GMatrixROmem5 extends DenseMatrix{self.M==self.N} {
    // Timer & PAPI performance counters
    val TIMER_TOTAL = 0;
    val TIMER_JMATRIX = 1;
    val TIMER_KMATRIX = 2;
    val TIMER_GENCLASS = 3;
    public val timer = new StatisticalTimer(4);

    // @Ifdef("__PAPI__") // XTENLANG-3132
    transient var papi:PAPI = new PAPI(); 
    val nPlaces:Long;

    // Standard conventional stuff
    private val bfs:BasisFunctions, mol:Molecule[QMAtom];
    val nOrbital:Long, norm:Rail[Double]; 
    val shellPairs:PlaceLocalHandle[Rail[ShellPair]];

    // RO stuff 
    val auxIntMat4K:DistDenseMatrix, halfAuxMat:DistDenseMatrix, distJ:DistDenseMatrix, distK:DistDenseMatrix;
    val auxIntMat4J:PlaceLocalHandle[Rail[Rail[Double]]], ylms:PlaceLocalHandle[Rail[Rail[Double]]];
    val ttemp4J:WorkerLocalHandle[Rail[Double]], ttemp4K:WorkerLocalHandle[Rail[Double]], taux:WorkerLocalHandle[Integral_Pack], tdk:WorkerLocalHandle[Rail[Double]];

    var roN:Int, roNK:Int, roL:Int; // 'var' because it can be overridden
    val roK:Long, roZ:Double, omega:Double, roThresh:Double;
    val funcAtPlace:Rail[Long], offsetAtPlace:Rail[Long];
    val dk:PlaceLocalHandle[Rail[Double]], e:PlaceLocalHandle[Rail[Double]];

    public def this(N:Long, bfs:BasisFunctions, mol:Molecule[QMAtom], nOrbital:Long, omega:Double,roThresh:Double):GMatrixROmem5{self.M==N,self.N==N} {     
        super(N, N);
        val jd = JobDefaults.getInstance();
        val nPlaces=Place.MAX_PLACES;

        Console.OUT.printf("\nGMatrixROmem5.x10 'public def this' %s...\n", getDateString());
        Console.OUT.printf("Omega=%5.3f thresh=%e rad=%7.3f\n", omega, roThresh, jd.rad);


        // Set up RO N, L, K, Z

        val l_n = new Rail[Int](jd.roN+3);
        val aux = new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ);
        if (omega>0.) { // long-range Ewald operator
            aux.getNL(l_n);
            roN=roNK=l_n(0);
            roL=l_n(roN+2);             
        } else { // full Coulomb operator
            roN=jd.roN;
            roL=jd.roL;
            if (jd.roNK==-1n) roNK=roN; else roNK=jd.roNK; 
        }
        val roK = (roL+1)*(roL+1);
        Console.OUT.printf("roN=%d roNK=%d roL=%d\n", roN, roNK, roL);

        // Some other handy variables
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;

        // Data distribution on a shell basis
        // Input: nPlaces, mol
        // Output: place2atom, place2func

        val npp=N/nPlaces;
        val place2atom = new Rail[Long](nPlaces+1);
        val place2func = new Rail[Long](nPlaces+1);
        val funcAtPlace = new Rail[Long](nPlaces);
        val offsetAtPlace = new Rail[Long](nPlaces+1);
        var placeID:Long = nPlaces-1; var func:Long=0;
        val noOfAtoms=mol.getNumberOfAtoms();

        place2atom(nPlaces)=noOfAtoms-1;
        place2func(nPlaces)=mol.getAtom(noOfAtoms-1).getBasisFunctions().size();

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
                     placeID--;
                }
            }
        }
        if (placeID>0n) {Console.ERR.println("Too many palaces!\n"); System.setExitCode(1n); throw new UnsupportedOperationException("Too many places! Last placeID: "+placeID);}      
        place2atom(0)=0; place2func(0)=0;
        for (ii in (0..nPlaces)) Console.OUT.printf("place %3d: atom=%5d function=%3d\n",ii,place2atom(ii),place2func(ii));

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

        var mu:Long=0,nu:Long=0,ind:Long=0,totFunc:Long=0,skip:Long=0;
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
                        for (var ii:Int=0n; ii<dConA; ii++) for (var jj:Int=0n; jj<dConB; jj++) 
                            contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2)
                                     /Math.pow(jd.roZ,aang+bang);  // See Szabo Ostlund 3.284-3.286 
                        contrib=Math.abs(contrib);
                        if (offsetAtPlace(placeID-1)<=nu && nu<offsetAtPlace(placeID) && mu > nu) skip++;
                        else if (contrib>=threshold) {
                            val maxL = new Rail[Int](roN+1); for (var ron:Long=0; ron<=roN; ron++) maxL(ron)=roL;
                            rawShellPairs(ind) = new ShellPair(aang,bang,aPoint,bPoint,zetaA,zetaB,conA,conB,dConA,dConB,mu,nu,maxL,contrib);     
                            ind++;
                            totFunc+=maxbraa*maxbrab; 
                        }  
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab; else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        val numSigShellPairs=place2ShellPair(placeID)=ind; 
        funcAtPlace(placeID-1) = mu-offsetAtPlace(placeID-1); 
        Console.OUT.printf("#basis functions1 = %5d #basis functions2 = %5d\n", funcAtPlace(placeID-1), totFunc);
        Console.OUT.printf("place %3d: offset = %5d shellpair #%5d\n", nPlaces,mu+1, place2ShellPair(nPlaces));
      
        Console.OUT.printf("Found %d significant shellpairs (skipped %d shellpairs)\n", numSigShellPairs+skip,skip);       

        /*for (ii in (0..nPlaces)) {
            Console.OUT.printf("place %3d: offset = %5d shellpair #%5d\n",ii,offsetAtPlace(ii),place2ShellPair(ii));
            Console.OUT.printf("#basis functions1 = %5d\n",funcAtPlace(ii),totFunc); 
        }*/


        val taux = new WorkerLocalHandle[Integral_Pack](() => new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        val shellPairs = PlaceLocalHandle.make[Rail[ShellPair]](
            PlaceGroup.WORLD,
            ()=>new Rail[ShellPair](place2ShellPair(here.id+1)-place2ShellPair(here.id),
                (i:Long) => rawShellPairs(i+place2ShellPair(here.id))));

        val roL_val = roL;
        val ylms = PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>new Rail[Rail[Double]](place2ShellPair(here.id+1)-place2ShellPair(here.id),
                (i:Long) => 
                {
                    val sp = shellPairs()(i), tempY = new Rail[Double](sp.dconA*sp.dconB*(sp.maxL+1)*(sp.maxL+1));
                    taux().genClassY(sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.dconA, sp.dconB, roL_val, tempY);
                    tempY
                } ));

        val auxIntMat4J = PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>new Rail[Rail[Double]](place2ShellPair(here.id+1)-place2ShellPair(here.id),
                (i:Long) => 
                {
                    val pid = here.id, sp=shellPairs()(i);
                    val mult=(Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                    val colStart=offsetAtPlace(pid);
                    val colStop=offsetAtPlace((pid+mult)%nPlaces);

                    val nu=sp.nu; 
                    var size:Long=0L;
                    if ( ( (colStart<colStop) && ((colStart<=nu) && (nu<colStop)) ) ||
                         ( (colStart>=colStop) && ( (nu<colStop) || (nu>=colStart) ) ) )
                         size=sp.maxbraa*sp.maxbrab*roK;
                    if (offsetAtPlace(pid)<=nu && nu<offsetAtPlace(pid+1) && sp.mu > sp.nu) size=0L;
                    new Rail[Double](size)
                } ));        

        val cbs_auxInt= new Rail[Long](1);  cbs_auxInt(0)=N*roK;
        val auxIntGrid4K = new Grid(funcAtPlace, cbs_auxInt);
        val auxIntMat4K = DistDenseMatrix.make(auxIntGrid4K);

        val cbs_HalfAuxInt= new Rail[Long](1);  cbs_HalfAuxInt(0)=nOrbital*roK;
        val halfAuxGrid = new Grid(funcAtPlace, cbs_HalfAuxInt);
        val halfAuxMat = DistDenseMatrix.make(halfAuxGrid);

        val ttemp4J = new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](maxam1*maxam1*roK));
        val ttemp4K= new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](maxam1*maxam1*roK));
        val tdk = new WorkerLocalHandle[Rail[Double]](() => new Rail[Double](roK));

        val cbs_nSquareMat = new Rail[Long](1); cbs_nSquareMat(0) = N;
        val nSquareMatGrid = new Grid(funcAtPlace, cbs_nSquareMat);
        val distJ = DistDenseMatrix.make(nSquareMatGrid);
        val distK = DistDenseMatrix.make(nSquareMatGrid);

        this.dk = PlaceLocalHandle.make[Rail[Double]](
            PlaceGroup.WORLD, 
            () => new Rail[Double](roK)
        ); 

        this.e = PlaceLocalHandle.make[Rail[Double]](
            PlaceGroup.WORLD, 
            () => new Rail[Double](2)
        );  
        
        this.norm = bfs.getNormalizationFactors(); // Vector
        this.roNK=roNK; this.roN=roN; this.roL=roL; this.roK = roK; this.roZ=jd.roZ;  this.nPlaces=nPlaces;
        this.distJ=distJ; this.distK=distK;        
        this.ttemp4J=ttemp4J; this.ttemp4K=ttemp4K; this.tdk=tdk; this.taux = taux; 
        this.auxIntMat4J=auxIntMat4J; this.auxIntMat4K=auxIntMat4K;  this.halfAuxMat=halfAuxMat;
        this.offsetAtPlace=offsetAtPlace; this.funcAtPlace=funcAtPlace;        
        this.bfs = bfs; this.mol = mol; this.nOrbital = nOrbital; this.omega=omega; this.roThresh=roThresh;
        this.ylms = ylms; this.shellPairs = shellPairs;

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
        val timer=this.timer, nPlaces=this.nPlaces; 
        val shellPairs=this.shellPairs; val ylms=this.ylms;
        val N=this.N, nOrbital=this.nOrbital, roN=this.roN, roNK=this.roNK, roL=this.roL, roK=this.roK, roZ=this.roZ, norm=this.norm, e = this.e;
        val funcAtPlace=this.funcAtPlace, offsetAtPlace=this.offsetAtPlace, taux=this.taux, dk=this.dk;
        val auxIntMat4K=this.auxIntMat4K, halfAuxMat=this.halfAuxMat, auxIntMat4J=this.auxIntMat4J, ttemp4J=this.ttemp4J, ttemp4K=this.ttemp4K, tdk=this.tdk;
        val distJ=this.distJ, distK=this.distK; 
        timer.start(TIMER_TOTAL); 
        val jd = JobDefaults.getInstance();

        finish ateach(place in Dist.makeUnique()) {
            val pid = here.id;
            // Console.OUT.println("pid=" + pid + " starts...");

            // For faster access 
            val shp=shellPairs(), ylmp = ylms();
            val localAuxJ=auxIntMat4J(), localMatK=auxIntMat4K.local(); 
            val localJ=distJ.local(), localK=distK.local();
            val dkp = dk(), ep=e(); 
            
            ep.clear(); localJ.reset(); localK.reset();

            for (ron in 0n..roN) {
                // Console.OUT.println("Aux - distributed ron="+ron);
                timer.start(TIMER_GENCLASS); 
                dkp.clear();
                tdk.applyLocal((d:Rail[Double]) => { d.clear(); });

                // Aux & D
                finish DivideAndConquerLoop1D(0, shp.size).execute(
                (spInd:Long) => {
                    val sp=shp(spInd);
                    val maxLron=sp.L(ron);
                    if (maxLron>=0) {
                        val aux = taux();
                        var tempJ:Rail[Double]=localAuxJ(spInd);
                        if (tempJ.size==0L) tempJ=ttemp4J();                     

                        val tempK=ttemp4K(); val myThreaddk=tdk();
                        val maxLm=(maxLron+1)*(maxLron+1);
                        val y=ylmp(spInd);
                        aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, tempJ, ron, maxLron, y, sp.maxL);

                        var ind:Long=0;
                        val musize=sp.mu2-sp.mu+1; val nusize=sp.nu2-sp.nu+1;

                        if (localAuxJ(spInd).size==0L) {
                            for (var tmu:Long=sp.mu,ttmu:Long=0; tmu<=sp.mu2; tmu++,ttmu++) for (var tnu:Long=sp.nu,ttnu:Long=0; tnu<=sp.nu2; tnu++,ttnu++) {
                                val nrm=norm(tmu)*norm(tnu);
                                for (var rolm:Long=0; rolm<maxLm; rolm++) {
                                    val normAux = (tempJ(ind++)*=nrm);       
                                    tempK((ttnu*roK+rolm)*musize+ttmu) = normAux;  
                                }
                            }
                        } else if (sp.mu!=sp.nu) {
                            for (var tmu:Long=sp.mu,ttmu:Long=0; tmu<=sp.mu2; tmu++,ttmu++) for (var tnu:Long=sp.nu,ttnu:Long=0; tnu<=sp.nu2; tnu++,ttnu++) {
                                val scdmn=density(tmu,tnu); val nrm=norm(tmu)*norm(tnu); 
                                for (var rolm:Long=0; rolm<maxLm; rolm++) {
                                    val normAux = (tempJ(ind++)*=nrm);       
                                    myThreaddk(rolm) += 2.*scdmn*normAux; 
                                    tempK((ttnu*roK+rolm)*musize+ttmu) = normAux;  
                                }
                            }
                        } else {
                            for (var tmu:Long=sp.mu,ttmu:Long=0; tmu<=sp.mu2; tmu++,ttmu++) for (var tnu:Long=sp.nu,ttnu:Long=0; tnu<=sp.nu2; tnu++,ttnu++) {
                                val scdmn=density(tmu,tnu); val nrm=norm(tmu)*norm(tnu); 
                                for (var rolm:Long=0; rolm<maxLm; rolm++) {
                                    val normAux = (tempJ(ind++)*=nrm);       
                                    myThreaddk(rolm) += scdmn*normAux; 
                                    tempK((ttnu*roK+rolm)*musize+ttmu) = normAux;  
                                }
                            }
                        }
                        ind=0;
                        val rows = sp.mu2-sp.mu+1;
                        for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) {
                            DenseMatrix.copySubset(tempK, ind, localMatK, sp.mu-offsetAtPlace(pid), tnu*roK, rows, maxLm);
                            ind += rows*maxLm;
                        }

                        if (offsetAtPlace(pid)<=sp.nu && sp.nu<offsetAtPlace(pid+1) && sp.mu < sp.nu) {
                            ind=0;
                            for (var tmu:Long=sp.mu,ttmu:Long=0; tmu<=sp.mu2; tmu++,ttmu++) for (var tnu:Long=sp.nu,ttnu:Long=0; tnu<=sp.nu2; tnu++,ttnu++) for (var rolm:Long=0; rolm<maxLm; rolm++) 
                                tempK((ttmu*roK+rolm)*nusize+ttnu) = tempJ(ind++); 
                            val rows2 = sp.nu2-sp.nu+1; ind=0;
                            for (var tmu:Long=sp.mu; tmu<=sp.mu2; tmu++) {
                                DenseMatrix.copySubset(tempK, ind, localMatK, sp.nu-offsetAtPlace(pid), tmu*roK, rows2, maxLm);
                                ind += rows2*maxLm;
                            }
                        }                                    
                            

                    }
                }
                );

                tdk.reduceLocal(dkp, (a:Rail[Double],b:Rail[Double]) => RailUtils.map(a, b, a, (x:Double,y:Double)=>x+y));
                timer.stop(TIMER_GENCLASS);

                finish {
                    async Team.WORLD.allreduce[Double](dkp, 0L, dkp, 0L, dkp.size, Team.ADD);

                    timer.start(TIMER_KMATRIX);
                    if (ron <= roNK) {
                        val A=new DenseMatrix(funcAtPlace(pid)*roK, N, auxIntMat4K.local().d);
                        val B=new DenseMatrix(funcAtPlace(pid)*roK, nOrbital, halfAuxMat.local().d); 
                        DenseMatrixBLAS.compMultTrans(A, mos, B, [funcAtPlace(pid)*roK, nOrbital, N], false);
                        //cannot do B.multTrans(A, mos, false); -  mos is [N,N] rather than [nObital,N]
                    }
                    timer.stop(TIMER_KMATRIX);
                }

                if (ron <= roNK) { // This produces K/2
                    timer.start(TIMER_KMATRIX);
                    // Console.OUT.println("K - distributed"); 
                    val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                    val a=halfAuxMat.local();                   
                    for (var blk:Long=0; blk<mult; blk++) {
                        val qid=(pid+blk)%nPlaces;
                        val b=at(Place(qid)) {halfAuxMat.local()};
                        val noff=offsetAtPlace(qid);
                        val c=new DenseMatrix(a.M,b.M);
                        c.multTrans(a, b, false);
                        for (var j:Long=0; j<c.N; j++) for (var i:Long=0; i<c.M; i++) {
                            localK(i,noff+j)+=c(i,j);
                        }
                    }
                    timer.stop(TIMER_KMATRIX);
                }

                timer.start(TIMER_JMATRIX); 
                // Console.OUT.println("J - distributed"); 
                finish DivideAndConquerLoop1D(0, shp.size).execute(
                (spInd:Long) => {
                    val sp=shp(spInd);                
                    val temp=localAuxJ(spInd);
                    val maxLron = sp.L(ron);
                    if (temp.size>0 && maxLron>=0n) {
                        val maxLm=(maxLron+1)*(maxLron+1); var ind:Long=0L; 
 
                        /*if (sp.mu!=sp.nu && offsetAtPlace(pid)<=sp.nu && sp.nu<offsetAtPlace(pid+1)) {
                            for (var tmu:Long=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) {
                                var jContrib:Double=0.;  
                                val tmuoff = tmu-offsetAtPlace(pid);
                                val tnuoff = tnu-offsetAtPlace(pid);
                                for (var rolm:Long=0; rolm<maxLm; rolm++) 
                                    jContrib += dkp(rolm)*temp(ind++);
                                localJ(tnuoff,tmu)+= jContrib;
                                localJ(tmuoff,tnu)+= jContrib;
                            } 
                        } else {*/
                            for (var tmu:Long=sp.mu, tmuoff:Long=tmu-offsetAtPlace(pid); tmu<=sp.mu2; tmu++,tmuoff++) {
                                for (var tnu:Long=sp.nu; tnu<=sp.nu2; tnu++) {
                                    var jContrib:Double=0.;
                                    for (var rolm:Long=0; rolm<maxLm; rolm++) 
                                        jContrib += dkp(rolm)*temp(ind++);
                                    localJ(tmuoff,tnu) += jContrib;
                                }
                            }
                        //}
                    }
                }
                );
                timer.stop(TIMER_JMATRIX);

                Team.WORLD.barrier();
            }

            //Console.OUT.printf("\nG matrix\n");

            // Fix J
            finish for (spInd in 0..(shp.size-1)) async {
                val sp=shp(spInd);                
                if (localAuxJ(spInd).size>0L && offsetAtPlace(pid)<=sp.nu && sp.nu<offsetAtPlace(pid+1) && sp.mu!=sp.nu) {                     
                    for (var tmu:Long=sp.mu,tmuoff:Long = tmu-offsetAtPlace(pid); tmu<=sp.mu2; tmu++,tmuoff++) {
                        for (var tnu:Long=sp.nu, tnuoff:Long = tnu-offsetAtPlace(pid); tnu<=sp.nu2; tnu++,tnuoff++) {
                            localJ(tnuoff,tmu) = localJ(tmuoff,tnu);
                        }
                    } 
                }
            }

            // These variables are used in the two sections below:
            val rowCount=localJ.M; val colCount=localJ.N;
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace(pid);
            val colStop=offsetAtPlace((pid+mult)%nPlaces);

            // Calculate eJ and eK
            for (var j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j++) for (var i:Long=0,ii:Long=offsetAtPlace(pid); i<rowCount; i++,ii++)
                    {ep(0)-=.5*density(ii,j)*localJ(i,j); ep(1)-=.5*density(ii,j)*localK(i,j);}
            if (colStart<colStop) {
                for (var j:Long=colStart; j<colStop; j++) for (var i:Long=0,ii:Long=offsetAtPlace(pid); i<rowCount; i++,ii++)
                    {ep(0)+=density(ii,j)*localJ(i,j); ep(1)+=density(ii,j)*localK(i,j);}
            } else {
                for (var j:Long=colStart; j<colCount; j++) for (var i:Long=0,ii:Long=offsetAtPlace(pid); i<rowCount; i++,ii++)
                    {ep(0)+=density(ii,j)*localJ(i,j); ep(1)+=density(ii,j)*localK(i,j);}
                for (var j:Long=0; j<colStop; j++) for (var i:Long=0,ii:Long=offsetAtPlace(pid); i<rowCount; i++,ii++)
                    {ep(0)+=density(ii,j)*localJ(i,j); ep(1)+=density(ii,j)*localK(i,j);}
            }
            Team.WORLD.allreduce[Double](ep, 0L, ep, 0L, ep.size, Team.ADD);
            if (here == Place.FIRST_PLACE) Console.OUT.printf("  EJ = %.10f EK=%.10f\n", ep(0)/jd.roZ, -0.5*ep(1)/jd.roZ);

            // Combine J and K to form G (stored in J)
            if (colStart<colStop) {
                for (var j:Long=colStart; j<colStop; j++) for (var i:Long=0; i<rowCount; i++)
                    localJ(i,j)-=localK(i,j);
            } else {
                for (var j:Long=colStart; j<colCount; j++) for (var i:Long=0; i<rowCount; i++)
                    localJ(i,j)-=localK(i,j);
                for (var j:Long=0; j<colStop; j++) for (var i:Long=0; i<rowCount; i++)
                    localJ(i,j)-=localK(i,j);
            }

            // Report time
            Team.WORLD.allreduce[Long](timer.total, 0L, timer.total, 0L, timer.total.size, Team.MAX);
            if (here == Place.FIRST_PLACE) {
                val tINT = (timer.total(TIMER_GENCLASS) as Double)/1e9;
                val tJ = (timer.total(TIMER_JMATRIX) as Double)/1e9;
                val tK = (timer.total(TIMER_KMATRIX) as Double)/1e9;
                Console.OUT.printf("Time INT = %.2f s J = %.2f s K = %.2f s\n", tINT, tJ, tK);
                Console.OUT.flush();
            }
        }

        // Copy G to place 0 // can improve further by copying only contributing blocks
        for (pid in 0..(nPlaces-1)) {
            val mat=at(Place(pid)) {distJ.local()};
            DenseMatrix.copySubset(mat, 0, 0, this, offsetAtPlace(pid), 0, mat.M, mat.N);
        }

        // Fix G        
        for (var pid:Long=0; pid<nPlaces; pid++) {
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            for (var qid:Long=mult+pid; qid<nPlaces+pid; qid++) {
                val qq=qid%nPlaces;
                for (var i:Long=offsetAtPlace(pid); i<offsetAtPlace(pid+1); i++) for (var j:Long=offsetAtPlace(qq); j<offsetAtPlace(qq+1); j++)                
                    this(i,j) = this(j,i); 
            }
        }

        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);
 
        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }      
    }

    private static struct DivideAndConquerLoop1D(start:Long, end:Long) {
        // TODO allow grain size > 1
        public def this(start:Long, end:Long) {
            property(start, end);
        }

        public def execute(body:(idx:Long) => void) {
            if ((end-start) > 1L) {
                val firstHalf = DivideAndConquerLoop1D(start, (start+end)/2L);
                val secondHalf = DivideAndConquerLoop1D((start+end)/2L, end);
                async firstHalf.execute(body);
                secondHalf.execute(body);
            } else {
                body(start);
            }
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
