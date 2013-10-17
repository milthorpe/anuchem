/*  This file is part of ANUChem.
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
import x10.compiler.NonEscaping;
import x10.io.IOException;
import x10.lang.PlaceLocalHandle;
import x10.matrix.blas.DenseMatrixBLAS;
import x10.matrix.block.Grid;
import x10.matrix.DenseMatrix;
import x10.matrix.dist.DistDenseMatrix;
import x10.regionarray.Dist;
import x10.util.GrowableRail;
import x10.util.Pair;
import x10.util.RailUtils;
import x10.util.Team;
import x10.util.WorkerLocalHandle;

import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ro.Integral_Pack;
import au.edu.anu.qm.ShellPair; 
import au.edu.anu.util.StatisticalTimer;

import edu.utk.cs.papi.PAPI;

@NativeCPPInclude("mkl_math.h")
@NativeCPPInclude("omp.h")

public class GMatrixROmem5 extends DenseMatrix{self.M==self.N} {
    // Timer & PAPI performance counters
    public val timer=new StatisticalTimer(7);
    val TIMER_TOTAL=0;
    val TIMER_AUX=1;
    val TIMER_JMATRIX=2;
    val TIMER_K=3;
    val TIMER_DSYRK=4;
    val TIMER_DGEMM=5;
    val TIMER_COP=6;

    transient var papi:PAPI=new PAPI(); // @Ifdef("__PAPI__") // XTENLANG-3132

    private val bfs:BasisFunctions, mol:Molecule[QMAtom];

    val halfAuxMat:DistDenseMatrix, distJ:DistDenseMatrix, distK:DistDenseMatrix;
    val auxIntMat4J:PlaceLocalHandle[Rail[Rail[Double]]], auxIntMat4K:PlaceLocalHandle[Rail[Long]], remoteBlockK:PlaceLocalHandle[RemoteBlock], ylms:PlaceLocalHandle[Rail[Rail[Double]]], shellPairs:PlaceLocalHandle[Rail[ShellPair]];
    val dk:PlaceLocalHandle[Rail[Double]], e:PlaceLocalHandle[Rail[Double]], shellPairRange:PlaceLocalHandle[Rail[Long]]; 
    val ttemp4K:WorkerLocalHandle[Rail[Double]], taux:WorkerLocalHandle[Integral_Pack], tdk:WorkerLocalHandle[Rail[Double]], tB1:WorkerLocalHandle[DenseMatrix];
    val nOrbitals:Long, norm:Rail[Double], roN:Int, roNK:Int, roL:Int, roK:Int; 
    val roZ:Double, omega:Double, roThresh:Double, shellAtPlace:Rail[Long], funcAtPlace:Rail[Long], offsetAtPlace:Rail[Long];

    public def this(N:Long, bfs:BasisFunctions, mol:Molecule[QMAtom], nOrbitals:Long, omega:Double, roThresh:Double):GMatrixROmem5{self.M==N, self.N==N} {     
        super(N, N);
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def this' %s...\n", getDateString());
        val timer=new StatisticalTimer(1), jd=JobDefaults.getInstance(), nPlaces=Place.MAX_PLACES, nAtoms=mol.getNumberOfAtoms(), 
            maxam=bfs.getShellList().getMaximumAngularMomentum(), maxam1=(maxam+1)*(maxam+2)/2,
            l_n=new Rail[Int](jd.roN+3), aux=new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ),
            shellAtPlace=new Rail[Long](nPlaces), funcAtPlace=new Rail[Long](nPlaces), offsetAtPlace=new Rail[Long](nPlaces+1), place2atom=new Rail[Long](nPlaces+1), place2func=new Rail[Long](nPlaces+1);
        var nShells:Long=0, mu:Long=0, nu:Long=0, ind:Long=0, totY:Long=0, totJ:Long=0, skip:Long=0, pid:Long=nPlaces-1, func:Long=0;

        timer.start(0);
        // Set up nShells and RO variables for later use
        for (var a:Long=0; a<nAtoms; nShells+=mol.getAtom(a++).getBasisFunctions().size()); 
        if (omega>0.) { // long-range Ewald operator
            aux.getNL(l_n);
            roN=l_n(0);
            roL=l_n(roN+2);             
        } else { // full Coulomb operator
            roN=jd.roN;
            roL=jd.roL;
        }
        if (jd.roNK==-1n || jd.roNK>roN) roNK=roN; else roNK=jd.roNK;  // default/not provided roNK
        val roK=(roL+1n)*(roL+1n);

        // Preliminary work sharing: divide first by atoms, then by basis
        // functions so that each place has approximately equal total
        // angular momentum over all functions. Run backward so that 
        // (if there is load imbalance) the head node has less work.
        place2atom(nPlaces)=nAtoms-1; place2atom(0)=0;
        place2func(nPlaces)=mol.getAtom(nAtoms-1).getBasisFunctions().size(); place2func(0)=0;
        val npp=N/nPlaces; 
        for(var a:Long=nAtoms-1; a>=0 && pid>0; a--) { // centre a  
            val aFunc=mol.getAtom(a).getBasisFunctions(), naFunc=aFunc.size();          
            for(var i:Long=naFunc-1; i>=0 && pid>0; i--) { // basis functions on a
                val iaFunc=aFunc.get(i), aa=iaFunc.getTotalAngularMomentum();
                func+=(aa+1)*(aa+2)/2;
                if (N-func<pid*npp) {
                     place2atom(pid)=a;
                     place2func(pid)=i;
                     pid--;
                }
            }
        }
        if (pid>0n) {
            Console.ERR.println("Too many places! Last pid: "+pid); System.setExitCode(1n); 
            throw new UnsupportedOperationException("Too many places! Last pid: "+pid);
        }      
        pid=1; mu=0; // Run forward
        offsetAtPlace(0)=0; offsetAtPlace(nPlaces)=N; 
        for(var a:Long=0; a<nAtoms; a++) { // centre a  
            val aFunc=mol.getAtom(a).getBasisFunctions(), naFunc=aFunc.size();          
            for(var i:Long=0; i<naFunc; i++) { // basis functions on a
                val iaFunc=aFunc.get(i), aa=iaFunc.getTotalAngularMomentum();                   
                if (a==place2atom(pid) && i==place2func(pid)) {
                     funcAtPlace(pid-1)=mu-offsetAtPlace(pid-1);
                     offsetAtPlace(pid)=mu;
                     pid++;
                }
                shellAtPlace(pid-1)++;
                mu+=(aa+1)*(aa+2)/2;
            }
        }
        funcAtPlace(pid-1)=mu-offsetAtPlace(pid-1); 

        @Ifdef("__DEBUG__") {
            Console.OUT.printf("roN=%d roNK=%d roL=%d Omega=%.3f thresh=%e rad=%.3f\n", roN, roNK, roL, omega, roThresh, jd.rad);
            Console.OUT.printf("nAtoms=%d nShells=%d N=%d maxam=%d maxam1=%d\n\n", nAtoms, nShells, N, maxam, maxam1);
            Console.OUT.println("Preliminary work division by atom and function shell");
            for (i in (0..(nPlaces-1))) Console.OUT.printf("Place %3d: Atom=%5d Function=%3d #Shell=%d\n", i, place2atom(i), place2func(i), shellAtPlace(i));
            timer.stop(0);
            Console.OUT.println ("    GMatrixROmem5 Initialization 'Initial Assessment' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
            Console.OUT.printf("\n");
            timer.start(0);
        }

        // distributed generation of shellPairs
        val threshold=roThresh*jd.roZ*jd.roZ*1e-3; 
        // ** Threshold must be relative to roThresh *** otherwise Z scaling will cause a problem: This is effectively a density threshold RO Thesis (2.26)    
        val roL_val=roL,roN_val=roN,roZ_val=jd.roZ,nShells_val=nShells, sizeInfo=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](4));
        val shellPairRange=PlaceLocalHandle.make[Rail[Long]](PlaceGroup.WORLD, ()=>new Rail[Long](shellAtPlace(here.id)));
        val shellPairs=PlaceLocalHandle.make[Rail[ShellPair]](
            PlaceGroup.WORLD, 
            ()=> {            
            val pid=here.id, info=sizeInfo(), range=shellPairRange();
            val localShellPairs = new GrowableRail[ShellPair](nShells_val*nShells_val);
            var mu:Long = offsetAtPlace(pid), nu:Long=0, shell:Long=0;
            for (a in place2atom(pid)..place2atom(pid+1)) { // centre a  
                val aFunc = mol.getAtom(a).getBasisFunctions();
                val naFunc = aFunc.size();
                // basis functions on a: careful/tricky
                val minFunc = (a==place2atom(pid)) ? place2func(pid) : 0;
                val maxFunc = (a==place2atom(pid+1)) ? place2func(pid+1)-1 : naFunc-1;
                for (i in minFunc..maxFunc) {
                    val iaFunc = aFunc.get(i);    
                    val aa = iaFunc.getTotalAngularMomentum();
                    val aang = iaFunc.getTotalAngularMomentum();
                    val aPoint = iaFunc.origin;
                    val zetaA = iaFunc.exponents;
                    val conA = iaFunc.coefficients; 
                    val maxbraa = (aa+1)*(aa+2)/2;
                    finish for (b in 0..(nAtoms-1)) { // centre b
                        val bFunc = mol.getAtom(b).getBasisFunctions();
                        val nbFunc = bFunc.size();
                        for (j in 0..(nbFunc-1)) { // basis functions on b
                            val jbFunc = bFunc(j);
                            val bb = jbFunc.getTotalAngularMomentum();
                            val maxbrab = (bb+1)*(bb+2)/2; 
                            val spMu = mu;
                            val spNu = nu;
                            async {
                                val bang = jbFunc.getTotalAngularMomentum();
                                val bPoint = jbFunc.origin; 
                                val zetaB = jbFunc.exponents; 
                                val conB = jbFunc.coefficients; 
                                var contrib:Double = 0.; // conservative estimate from ss
                                val R2 = aPoint.distanceSquared(bPoint);
                                for (ii in 0..(conA.size-1)) {
                                    for (jj in 0..(conB.size-1)) {
                                        // See Szabo Ostlund 3.284-3.286
                                        contrib += conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2)/Math.pow(roZ_val, aang+bang);
                                    }
                                }
                                contrib=Math.abs(contrib); 
                                if (/*offsetAtPlace(pid) <= spNu && spNu < offsetAtPlace(pid+1) && spMu > spNu/*/false) {
                                    atomic info(3)++; 
                                } else if (contrib >= threshold) {
                                    val maxL = new Rail[Int](roN_val+1, roL_val); // change roL_val to smaller number
                                    val sp = new ShellPair(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, spMu, spNu, maxL, contrib);
                                    atomic {
                                        localShellPairs.add(sp);
                                        info(0)+=maxbraa*maxbrab; 
                                        info(1)+=conA.size*conB.size*roK;
                                        info(2)+=maxbraa*maxbrab*roK;
                                    }
                                } else atomic info(3)++;
                            }  
                            if (b!=nAtoms-1 || j!=nbFunc-1) nu += maxbrab;
                            else {mu += maxbraa; nu = 0;}
                        }    
                    }
                    range(shell++)=localShellPairs.size();
                    // Console.OUT.println(here + " Shell: " + (shell-1) + " ShellPair ID: " + range(shell-1));
                }   
            }    
            localShellPairs.toRail()
        });

        // Calculating various measure of load (im)balance
        var max:Double=funcAtPlace(0), min:Double=max, tot:Double=N, ideal:Double=tot/nPlaces, tot2:Double=0.;
        Console.OUT.printf("1. Number of mu at each place\nPlace  (Offset)  Functions  Fraction\n");
        for (i in (0..(nPlaces-1))) {
            val cost=funcAtPlace(i);                      
            max=Math.max(cost,max); min=Math.min(cost,min); tot2+=cost; 
            Console.OUT.printf("%5d  (%6d)  %9d  %7.2f%%\n", i, offsetAtPlace(i), cost, cost*100./tot);
        }
        Console.OUT.printf("Fractions add up to %.2f %%\n", tot2/tot*100.);
        Console.OUT.printf("Fraction of N at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %%\n\n", ideal, max, min, (max/ideal-1.)*100.);
        val maxRow = max as Long; // This is used for remoteBlockK allocation later

        tot=N*N*(nPlaces+1.)*.5/nPlaces; tot2=0.; ideal=max=min=tot/nPlaces;
        Console.OUT.printf("2. Number of mu nu (block of G/J/K) at each place (bad for X10_NPLACES=even number)\nPlace  Functions  Fraction\n");      
        for (i in (0..(nPlaces-1))) {
            val mult=(Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && i<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace(i), colStop=offsetAtPlace((i+mult)%nPlaces);
            val col=colStart<colStop?colStop-colStart:colStop+N-colStart;
            val cost=funcAtPlace(i)*col;
            max=Math.max(cost,max); min=Math.min(cost,min); tot2+=cost;
            Console.OUT.printf("%5d  %9d  %7.2f%%\n", i, cost, cost*100./tot);            
        }
        Console.OUT.printf("Fractions add up to %.2f %% (due to rounding of N/nPlaces)\n",tot2/tot*100.);
        Console.OUT.printf("Block size at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %%\n\n", ideal, max, min, (max/ideal-1.)*100.);

        tot=N*(1.+N)*.5; tot2=0.; ideal=tot/nPlaces; max=min=at(Place(0)) sizeInfo()(0);
        Console.OUT.printf("3. Number of Aux(mu,nu) calculated at each place\nPlace  Functions  Fraction\n");    
        for (i in (0..(nPlaces-1))) {
            val cost=at(Place(i)) sizeInfo()(0);
            totY+=at(Place(i)) sizeInfo()(1);
            totJ+=at(Place(i)) sizeInfo()(2);
            skip+=at(Place(i)) sizeInfo()(3);
            max=Math.max(cost,max); min=Math.min(cost,min); tot2+=cost;
            Console.OUT.printf("%5d %10.0f  %7.2f%%\n", i, cost, cost*100./tot); 
        }
        Console.OUT.printf("Fractions add up to %.2f %% (due to rounding of N/nPlaces, granularity of shellpairs and shellpair cut-off)\n",tot2/tot*100.);
        Console.OUT.printf("Aux/D calculations at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %% (based on ideal)\n", ideal, max, min, (max/ideal-1.)*100.);
        ideal=tot2/nPlaces;
        Console.OUT.printf(" Imbalance cost=%.2f %% (adjusted), %d shellpairs skipped\n\n", (max/ideal-1.)*100., skip);

        Console.OUT.printf("Matrices size in MBs/64-bit double/\nJ, K, G, density, mos\t%.3f (each)\naux4J \t%.3f\nYlm  \t%.3f\naux4K \t%.3f\nhalfAux\t%.3f\n\n", N*N*8e-6, totJ*8e-6, totY*8e-6, N*N*roK*8e-6, nOrbitals*N*roK*8e-6);        

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to ShellPair' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);
        
        val taux=new WorkerLocalHandle[Integral_Pack](()=> new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        this.ylms=PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs();
                val ylms = new Rail[Rail[Double]](shp.size); 
                finish for (i in 0..(shp.size-1)) async {
                    val sp = shp(i);
                    val tempY = new Rail[Double](sp.conA.size * sp.conB.size * (sp.maxL+1)*(sp.maxL+1));
                    taux().genClassY(sp.bPoint, sp.aPoint, sp.zetaB, sp.zetaA, roL_val, tempY);
                    ylms(i) = tempY;
                } 
                ylms
            }
        );

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to ylms' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        this.auxIntMat4J=PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs(), pid = here.id;
                val auxIntMat4J = new Rail[Rail[Double]](shp.size);
                val mult = (Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                val colStart = offsetAtPlace(pid);
                val colStop = offsetAtPlace((pid+mult)%nPlaces);
                finish for (i in 0..(shp.size-1)) async {
                    val sp = shp(i), nu = sp.nu, mu = sp.mu;
                    var size:Long = 0;
                    if ( ( (colStart<colStop) && ((colStart<=nu) && (nu<colStop)) ) ||
                         ( (colStart>=colStop) && ( (nu<colStop) || (nu>=colStart) ) ) )
                        size = sp.maxbraa*sp.maxbrab*roK;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu > nu) size=0;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu < nu) size*=2;
                    auxIntMat4J(i) = new Rail[Double](size);
                }
                auxIntMat4J
            }
        );      

        this.auxIntMat4K=PlaceLocalHandle.make[Rail[Long]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs(), size=shp.size, pid=here.id;
                val auxIntMat4K = new Rail[Long](shp.size);
                finish for (i in 0..(size-1)) async {
                    val sp = shp(i), mu = sp.mu, nu = sp.nu;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu != nu) {
                        var j:Long=-1;
                        if (mu>nu) for (j=i-1; shp(j).mu!=nu || shp(j).nu!=mu; j--);
                        else for (j=i+1; shp(j).mu!=nu || shp(j).nu!=mu; j++);
                        auxIntMat4K(i)=j+size;
                        //if (shp(j).mu!=nu || shp(j).nu!=mu) Console.OUT.println ("i="+i+"j="+j);
                    } else  auxIntMat4K(i) = -1;
                }
                auxIntMat4K
            }
        );    

        val tbs = Math.max(maxRow*nOrbitals,maxam1*N)*roK; // if we use this for K too
        this.remoteBlockK=PlaceLocalHandle.make[RemoteBlock](PlaceGroup.WORLD, () => new RemoteBlock(tbs/*maxRow*nOrbitals*roK*/));

        val cbs_HalfAuxInt=new Rail[Long](1); cbs_HalfAuxInt(0)=nOrbitals*roK; val halfAuxIntGrid=new Grid(cbs_HalfAuxInt, funcAtPlace);
        this.halfAuxMat=DistDenseMatrix.make(halfAuxIntGrid);

        this.ttemp4K=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](maxam1*N*roK));
        this.tdk=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](roK));
        this.tB1=new WorkerLocalHandle[DenseMatrix](()=> new DenseMatrix(nOrbitals, roK*funcAtPlace(here.id)));

        val cbs_nSquareMat=new Rail[Long](1); cbs_nSquareMat(0)=N; val nSquareMatGrid=new Grid(funcAtPlace, cbs_nSquareMat);
        this.distJ=DistDenseMatrix.make(nSquareMatGrid); this.distK=DistDenseMatrix.make(nSquareMatGrid);

        this.dk=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](roK)); 
        this.e=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](2));
        this.shellPairRange=shellPairRange;

        this.shellAtPlace=shellAtPlace; this.offsetAtPlace=offsetAtPlace; this.funcAtPlace=funcAtPlace;
        this.roK=roK; this.roZ=jd.roZ; this.taux=taux; this.shellPairs=shellPairs;
        this.bfs=bfs; this.mol=mol; this.nOrbitals=nOrbitals; this.omega=omega; this.roThresh=roThresh; this.norm=bfs.getNormalizationFactors();

        timer.stop(0);
        Console.OUT.println("    GMatrixROmem5 Initialization 'total' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        @Ifdef("__PAPI__") { papi.initialize(); papi.countFlops(); papi.countMemoryOps(); }
        @Ifdef("__DEBUG__") { printShellInfo(); }
    }

    @Native("c++", "mkl_get_max_threads()") private native static def mklGetMaxThreads():Int;
    @Native("c++", "MKL_Set_Num_Threads(#a)") private native static def mklSetNumThreads(a:Int):void;

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def compute' %s...\n", getDateString()); 
        timer.start(TIMER_TOTAL); 

        val timer=this.timer; 
        val shellPairs=this.shellPairs; val ylms=this.ylms, distJ=this.distJ, distK=this.distK; 
        val N=this.N, nOrbitals=this.nOrbitals, roN=this.roN, roNK=this.roNK, roL=this.roL, roK=this.roK, roZ=this.roZ, norm=this.norm;
        val shellAtPlace=this.shellAtPlace, funcAtPlace=this.funcAtPlace, offsetAtPlace=this.offsetAtPlace, taux=this.taux;
        val dk=this.dk, e=this.e, shellPairRange=this.shellPairRange;
        val auxIntMat4K=this.auxIntMat4K, halfAuxMat=this.halfAuxMat, remoteBlockK=this.remoteBlockK;
        val auxIntMat4J=this.auxIntMat4J, ttemp4K=this.ttemp4K, tdk=this.tdk, tB1=this.tB1;
      
        finish ateach(place in Dist.makeUnique()) {
            val pid = here.id;
            val offsetHere = offsetAtPlace(pid);
            val nPlaces = Place.MAX_PLACES;

            val nThreads=Runtime.NTHREADS;
            GMatrixROmem5.setThread(nThreads);

            // For faster access
            val shp=shellPairs(), size=shp.size, ylmp=ylms();
            val localAuxJ=auxIntMat4J(), localAuxK=auxIntMat4K();
            val localJ=distJ.local(), localK=distK.local();
            val dkp=dk(), ep=e(), range=shellPairRange(); 

            val B1=new DenseMatrix(nOrbitals, roK*funcAtPlace(pid), halfAuxMat.local().d); 
            
            ep.clear(); localJ.reset(); localK.reset();
            
            for (ron in 0n..roN) {
                // Console.OUT.println("Aux - distributed ron="+ron);
                // Aux & D
                timer.start(TIMER_AUX);
                val doK = (ron <= roNK);
                dkp.clear();
                tdk.applyLocal((d:Rail[Double])=> { d.clear(); });
                for (i in 0..(shp.size-1)) if (0<=localAuxK(i) && localAuxK(i)<size) localAuxK(i)+=size;
                GMatrixROmem5.setThread(1n);
                
                finish DivideAndConquerLoop1D(0, shellAtPlace(pid)).execute(
                (sInd:Long)=> {
                    var spInd0:Long=0; 
                    if (sInd>0) spInd0=range(sInd-1);
                    val sp0=shp(spInd0);
                    val muSize=sp0.mu2-sp0.mu+1;
                    val tbk=ttemp4K(); 
                    tbk.clear();
                    val AuxMat = new DenseMatrix(roK*muSize, N, tbk);
                    val aux=taux(), myThreaddk=tdk();
                    for (var spInd:Long=spInd0; spInd<range(sInd); spInd++) {
                        val sp=shp(spInd);
                        val maxLron=sp.L(ron);
                        if (maxLron >= 0) {
                            val maxLm=(maxLron+1)*(maxLron+1), nuSize=sp.nu2-sp.nu+1, temp=localAuxJ(spInd), srcSpInd=localAuxK(spInd);
                            val off=(roK*muSize*sp.nu) as Int, asize=muSize*nuSize*roK;
                            var ind:Long=0; 
                            if (srcSpInd<0 || srcSpInd>=size) { // Generate Aux Ints and normalize ()
                                val y=ylmp(spInd);
                                aux.genClass(sp.bang, sp.aang, sp.bPoint, sp.aPoint, sp.zetaB, sp.zetaA, sp.conB, sp.conA, ron, maxLron, y, sp.maxL, off, tbk);
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val nrm=norm(mu)*norm(nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++) 
                                        tbk(off+ind)*=nrm;
                                }
                            } else if (doK || (srcSpInd>0 && sp.nu>sp.mu)) {  // Read Aux Ints (for K or J)
                                var src:Rail[Double] = localAuxJ(srcSpInd); // from remote location
                                if (src.size==0) { 
                                    src=localAuxJ(spInd); ind=asize;  // from local location
                                }
                                for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) 
                                    for (rolm in 0..(maxLm-1)) 
                                        tbk((nu*muSize+tmu)*roK+rolm)=src(ind++);
                            }

                            if (temp.size!=0) Rail.copy(tbk, off as Long, temp, 0, asize); // purely for J
                            else if (srcSpInd>=size) { //write to remote shp J - to be read and save time
                                val temp2=localAuxJ(srcSpInd-size); 
                                Rail.copy(tbk, off as Long, temp2, asize, asize);
                            } 
                            if (srcSpInd>=size) localAuxK(srcSpInd-size)-=size;

                            // J matter
                            ind=0;
                            if (sp.mu==sp.nu) { // for diagonal block of J
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val scdmn=density(mu, nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++)
                                        myThreaddk(rolm)+=scdmn*temp(ind); 
                                }
                            } else if (temp.size!=0) { // for the rest of J
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val scdmn=density(mu, nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++) 
                                        myThreaddk(rolm)+=2.*scdmn*temp(ind);
                                } 
                            }
                        }
                    }
                    if (doK) DenseMatrixBLAS.compMultTrans(mos, AuxMat, B1, [nOrbitals, AuxMat.M, N], [0, 0, 0, 0, 0, (sp0.mu-offsetHere)*roK], false);
                }
                );
                GMatrixROmem5.setThread(nThreads);
                tdk.reduceLocal(dkp, (a:Rail[Double], b:Rail[Double])=> RailUtils.map(a, b, a, (x:Double, y:Double)=>x+y));
                Team.WORLD.allreduce[Double](dkp, 0L, dkp, 0L, dkp.size, Team.ADD);
                timer.stop(TIMER_AUX);

                finish {
                    timer.start(TIMER_K);
                    if (doK) {
                        val remoteK = remoteBlockK();
                        val a = halfAuxMat.local();

                        val blocks = (Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                        var blk:Long = 1;
                        var nextBlockPlace:Long = -1;
                        finish {
                            if (blocks > 1) {
                                // overlap getting first remote block of K with DSYRK
                                nextBlockPlace = (pid+blk) % nPlaces;
                                async remoteK.fetchNext(halfAuxMat, nextBlockPlace);
                            }
                            timer.start(TIMER_DSYRK);
                            DenseMatrixBLAS.symRankKUpdateTrans(a, localK, [a.N, a.M], [0, 0, 0, offsetHere], true, true);
                            timer.stop(TIMER_DSYRK);
                        }

                        // This produces K/2 & TODO: ring broadcast ?
                        while (nextBlockPlace != -1) {
                            val thisBlockPlace = nextBlockPlace;
                            val thisBlock = remoteK.getCurrent();
                            blk++;
                            finish {
                                if (blk < blocks) {
                                    // overlap getting next remote block of K with DGEMM
                                    nextBlockPlace = (pid+blk) % nPlaces;
                                    async remoteK.fetchNext(halfAuxMat, nextBlockPlace);
                                } else {
                                    nextBlockPlace = -1;
                                }
                                timer.start(TIMER_DGEMM);
                                DenseMatrixBLAS.compTransMult(a, thisBlock, localK, [a.N, thisBlock.N, a.M], [0, 0, 0, 0, 0, offsetAtPlace(thisBlockPlace)], true);
                                timer.stop(TIMER_DGEMM);
                            }

                        }
                    }
                    timer.stop(TIMER_K);
                }

                timer.start(TIMER_JMATRIX); 
                // Console.OUT.println("J - distributed"); 
                finish DivideAndConquerLoop1D(0, shp.size, 16).execute(
                (spInd:Long)=> {
                    val sp=shp(spInd);
                    val auxJ=localAuxJ(spInd);
                    val maxLron=sp.L(ron);
                    if (auxJ.size > 0 && maxLron >=0n) {
                        val maxLm=(maxLron+1)*(maxLron+1);
                        var ind:Long=0; 
                        for (nu in sp.nu..sp.nu2) {
                            for (var mu:Long=sp.mu, muoff:Long=mu-offsetHere; mu<=sp.mu2; mu++, muoff++) {
                                var jContrib:Double=0.;
                                for (rolm in 0..(maxLm-1)) {
                                    jContrib += dkp(rolm) * auxJ(ind++);
                                }
                                localJ(muoff, nu) += jContrib;
                            }
                        }
                    }
                }
                );
                timer.stop(TIMER_JMATRIX);
                Team.WORLD.barrier();
            }

            //Console.OUT.printf("\nG matrix\n");
            // These variables are used in the two sections below:
            val rowCount=localJ.M;
            val colCount=localJ.N;
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace((pid+1)%nPlaces);
            val colStop=offsetAtPlace((pid+mult)%nPlaces);
            // Calculate eJ and eK
            for (var j0:Long=0, j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j0++, j++) {
                ep(0) += .5*density(j, j)*localJ(j0, j);
                ep(1) += .5*density(j, j)*localK(j0, j);
                for (var i0:Long=0, i:Long=offsetAtPlace(pid); i<j; i0++, i++) {
                    ep(0) += density(i, j)*localJ(i0, j);
                    ep(1) += density(i, j)*localK(i0, j);
                }
            }
            if (colStart < colStop) {
                for (var j:Long=colStart; j<colStop; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        ep(0) +=density(ii, j)*localJ(i, j);
                        ep(1) +=density(ii, j)*localK(i, j);
                    }
                }
            } else if (mult>1) {
                for (var j:Long=colStart; j<colCount; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        ep(0) +=density(ii, j)*localJ(i, j);
                        ep(1) +=density(ii, j)*localK(i, j);
                    }
                }
                for (var j:Long=0; j<colStop; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        ep(0) +=density(ii, j)*localJ(i, j);
                        ep(1) +=density(ii, j)*localK(i, j);
                    }
                }
            }
            Team.WORLD.allreduce[Double](ep, 0L, ep, 0L, ep.size, Team.ADD);
            if (here==Place.FIRST_PLACE) Console.OUT.printf("EJ= %.10f EK= %.10f\n", ep(0)/roZ, -0.5*ep(1)/roZ);

            // Combine J and K to form G (stored in J)
            for (var j0:Long=0, j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j0++, j++) 
                for (var i0:Long=0; i0<=j0; i0++) 
                    localJ(i0, j) -=localK(i0, j);
            if (colStart < colStop) {
                for (var j:Long=colStart; j<colStop; j++) {
                    for (var i:Long=0; i<rowCount; i++) {
                        localJ(i, j) -=localK(i, j);
                    }
                }
            } else if (mult>1) {
                for (var j:Long=colStart; j<colCount; j++) {
                    for (var i:Long=0; i<rowCount; i++) {
                        localJ(i, j) -=localK(i, j);
                    }
                }
                for (var j:Long=0; j<colStop; j++) {
                    for (var i:Long=0; i<rowCount; i++) {
                        localJ(i, j) -=localK(i, j);
                    }
                }
            }

            // Report time
            Team.WORLD.allreduce[Long](timer.total, 0L, timer.total, 0L, timer.total.size, Team.MAX);
            if (here==Place.FIRST_PLACE) {
                val tAux=(timer.total(TIMER_AUX) as Double)/1e9;
                val tJ=(timer.total(TIMER_JMATRIX) as Double)/1e9;
                val tDsyrk=(timer.total(TIMER_DSYRK) as Double)/1e9;
                val tDgemm=(timer.total(TIMER_DGEMM) as Double)/1e9;
                val tK=(timer.total(TIMER_K) as Double)/1e9;
                Console.OUT.printf("Time (seconds) Aux= %.2f J= %.2f K-DSYRK= %.2f K-DGEMM= %.2f K-TOT= %.2f\n", tAux, tJ, tDsyrk, tDgemm, tK);
                Console.OUT.flush();
            }
        }

        val nPlaces=Place.MAX_PLACES;
        // Copy G to place 0 // can improve further by copying only contributing blocks?
        for (pid in 0..(nPlaces-1)) {
            val mat=at(Place(pid)) { distJ.local() };
            DenseMatrix.copySubset(mat, 0, 0, this, offsetAtPlace(pid), 0, mat.M, mat.N);
        }
        // Fix G
        for (pid in 0..(nPlaces-1)) {
            for (var j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j++) 
                for (var i:Long=offsetAtPlace(pid); i<j; i++) 
                        this(j, i)=this(i, j);
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            for (var qid:Long=pid+1; qid<pid+mult; qid++) {
                val qq=qid%nPlaces;
                for (var i:Long=offsetAtPlace(pid); i<offsetAtPlace(pid+1); i++) {
                    for (var j:Long=offsetAtPlace(qq); j<offsetAtPlace(qq+1); j++) {
                        this(j, i)=this(i, j);
                    }
                }
            }
        }
        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9); 
        @Ifdef("__PAPI__"){ papi.printFlops(); papi.printMemoryOps();}
    }

    private static struct DivideAndConquerLoop1D(start:Long, end:Long, grainSize:Long) {
        public def this(start:Long, end:Long) {
            property(start, end, 1);
        }

        public def this(start:Long, end:Long, grainSize:Long) {
            property(start, end, grainSize);
        }

        public def execute(body:(idx:Long)=> void) {
            if ((end-start) > grainSize) {
                val firstHalf=DivideAndConquerLoop1D(start, (start+end)/2L, grainSize);
                val secondHalf=DivideAndConquerLoop1D((start+end)/2L, end, grainSize);
                async firstHalf.execute(body);
                secondHalf.execute(body);
            } else {
                for (i in start..(end-1)) {
                    body(i);
                }
            }
        }
    }

    /**
     * Adds matrix b to matrix a in parallel using all available threads.
     */
    private def parallelAdd(a:DenseMatrix, b:DenseMatrix){a.M==b.M,a.N==b.N} {
        val aRaw = a.d;
        val bRaw = b.d;
        val size = aRaw.size;
        val chunk = size / Runtime.NTHREADS;
        val remainder = size % Runtime.NTHREADS;
        finish for (t in 0..(Runtime.NTHREADS-1)) async {
            val start = (t < remainder) ? t*(chunk+1) : t*chunk + remainder;
            val end = ((t < remainder) ? (t+1)*(chunk+1) : (t+1)*chunk + remainder) - 1;
            for (i in start..end) {
                aRaw(i) += bRaw(i);
            }
        }
        return a;
    }

    /** 
     * This class supports a ring-like algorithm where a remote block of
     * a DistDenseMatrix is fetched and then some operation is performed
     * on it concurrently with the fetch of the next remote block.
     */
    private static class RemoteBlock {
        var currentData:Rail[Double];
        var nextData:Rail[Double];
        var currentDim:Pair[Long,Long];
        var nextDim:Pair[Long,Long];

        public def this(blockSize:Long) {
            currentData = new Rail[Double](blockSize);
            nextData = new Rail[Double](blockSize);
        }

        public def getCurrent():DenseMatrix {
            // swap next and current
            currentDim = nextDim;
            val temp = nextData;
            nextData = currentData;
            currentData = temp;
            return new DenseMatrix(currentDim.first, currentDim.second, currentData);
        }

        public def fetchNext(ddm:DistDenseMatrix, nextBlockPlace:Long) {
            // fetch next data from given place
            val nextDataRef = new GlobalRail(nextData);
            finish nextDim = at(Place(nextBlockPlace)) {
                val dataHere = ddm.local().d;
                Rail.asyncCopy(dataHere, 0, nextDataRef, 0, dataHere.size);
                Pair(ddm.local().M, ddm.local().N)
            };
        }
    }

    private def getDateString() {
        val result:String;
        try {
            val dateReader=Runtime.execForRead("date"); 
            result=dateReader.readLine();
        } catch (e:IOException) {
            // could not read date! use current time in milliseconds
            result=System.currentTimeMillis() + "ms";
        }
        return result;
    } 

    @Native("c++", "omp_get_num_threads()") private native static def ompGetNumThreads():Int;
    @Native("c++", "omp_set_num_threads(#a)") private native static def ompSetNumThreads(a:Int):void;

    private static def setThread(nT:Int) {
        @Ifdef("__MKL__") /*finish ateach(place in Dist.makeUnique())*/ { // not working?  better use -genv OMP_NUM_THREADS 4
            val t1=mklGetMaxThreads();
            val o1=ompGetNumThreads();
            ompSetNumThreads(nT);
            mklSetNumThreads(nT);
            val t2=mklGetMaxThreads();
            val o2=ompGetNumThreads();
            @Ifdef("__DEBUG__") { Console.OUT.println(here + ", mklGetMaxThreads() was " + t1 + " and is now set to " + t2 + " thread(s)."
                            + " ompGetNumThreads() was " + o1 + " and is now set to " + o2 + " thread(s)."); }   
        }
    }

    private @NonEscaping def printShellInfo() {
        finish ateach(place in Dist.makeUnique()) {
            val hostname=Runtime.execForRead("uname -n").readLine();      
            val np=Runtime.execForRead("echo $X10_NPLACES").readLine();
            val nt=Runtime.execForRead("echo $X10_NTHREADS").readLine();
            val gt=Runtime.execForRead("echo $GOTO_NUM_THREADS").readLine();
            val omp=Runtime.execForRead("echo $OMP_NUM_THREADS").readLine();
            Console.OUT.println(here + ", Runtime.NTHREADS=" + Runtime.NTHREADS + ", uname -n=" + hostname + ", X10_NPLACES="+np+", X10_NTHREADS="+nt+", GOTO_NUM_THREADS="+gt+ ", OMP_NUM_THREADS="+omp ); // if print out on separate lines, they can goes randomly.
            Console.OUT.flush();
        }
   }
}
