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

public class ROFockMethod(N:Long) {
    // Timer & PAPI performance counters
    public val timer=new StatisticalTimer(6);
    val TIMER_TOTAL=0;
    val TIMER_AUX=1;
    val TIMER_JMATRIX=2;
    val TIMER_K=3;
    val TIMER_DSYRK=4;
    val TIMER_DGEMM=5;

    transient var papi:PAPI=new PAPI(); // @Ifdef("__PAPI__") // XTENLANG-3132

    val halfAuxMat:DistDenseMatrix, distJ:DistDenseMatrix, distK:DistDenseMatrix;
    val auxJMat_plh:PlaceLocalHandle[Rail[Rail[Double]]], auxKMat_plh:PlaceLocalHandle[Rail[Long]], remoteBlockK_plh:PlaceLocalHandle[RemoteBlock], ylms_plh:PlaceLocalHandle[Rail[Rail[Double]]], shellPairs_plh:PlaceLocalHandle[Rail[ShellPair]];
    val dlm_plh:PlaceLocalHandle[Rail[Double]], e_plh:PlaceLocalHandle[Rail[Double]], shellPairRange_plh:PlaceLocalHandle[Rail[Long]]; 
    val auxK_wlh:WorkerLocalHandle[Rail[Double]], intPack_wlh:WorkerLocalHandle[Integral_Pack], dlm_wlh:WorkerLocalHandle[Rail[Double]];
    val nOrbitals:Long, norm:Rail[Double], roN:Int, roNK:Int, roK:Int; 
    val roZ:Double, shellAtPlace:Rail[Long], funcAtPlace:Rail[Long], offsetAtPlace:Rail[Long];

    public def this(N:Long, bfs:BasisFunctions, mol:Molecule[QMAtom], nOrbitals:Long, omega:Double, roThresh:Double) {     
        Console.OUT.printf("\nROFockMethod.x10 'public def this' %s...\n", getDateString());
        property(N);
        val maxam=bfs.getShellList().getMaximumAngularMomentum();
        val maxam1=(maxam+1)*(maxam+2)/2;
        val timer=new StatisticalTimer(1), jd=JobDefaults.getInstance(), nPlaces=Place.MAX_PLACES,
            shellAtPlace=new Rail[Long](nPlaces), funcAtPlace=new Rail[Long](nPlaces), offsetAtPlace=new Rail[Long](nPlaces+1);

        timer.start(0);
        // Set up nShells and RO variables for later use
        var nShells:Long = 0;
        for (atom in mol.getAtoms()) nShells += atom.getBasisFunctions().size(); 

        val roL:Int;
        if (omega>0.) { // long-range Ewald operator
            val aux = new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ);
            val l_n = new Rail[Int](jd.roN+3);
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
        val nAtoms=mol.getNumberOfAtoms();
        val place2atom=new Rail[Long](nPlaces+1);
        place2atom(nPlaces)=nAtoms-1; place2atom(0)=0;
        val place2func=new Rail[Long](nPlaces+1);
        place2func(nPlaces)=mol.getAtom(nAtoms-1).getBasisFunctions().size(); place2func(0)=0;
        val npp=N/nPlaces;
        var pid:Long=nPlaces-1;
        var func:Long=0;
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
        if (pid > 0) {
            Console.ERR.println("Too many places! Last pid: "+pid); System.setExitCode(1n); 
            throw new UnsupportedOperationException("Too many places! Last pid: "+pid);
        }

        var mu:Long=0;
        pid=1;// Run forward
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
            Console.OUT.println ("    ROFockMethod Initialization 'Initial Assessment' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
            Console.OUT.printf("\n");
            timer.start(0);
        }

        // distributed generation of shellPairs_plh
        val threshold=roThresh*jd.roZ*jd.roZ*1e-3; 
        // ** Threshold must be relative to roThresh *** otherwise Z scaling will cause a problem: This is effectively a density threshold RO Thesis (2.26)    
        val roN_val=roN,roZ_val=jd.roZ,nShells_val=nShells, sizeInfo=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](4));
        val shellPairRange_plh=PlaceLocalHandle.make[Rail[Long]](PlaceGroup.WORLD, ()=>new Rail[Long](shellAtPlace(here.id)));
        val shellPairs_plh=PlaceLocalHandle.make[Rail[ShellPair]](
            PlaceGroup.WORLD, 
            ()=> {            
            val pid=here.id, info=sizeInfo(), range=shellPairRange_plh();
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
                                if (contrib >= threshold) {
                                    val maxL = new Rail[Int](roN_val+1, roL); // change roL to smaller number
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
        val maxRow = max as Long; // This is used for remoteBlockK_plh allocation later

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
        var totY:Long=0, totJ:Long=0, skip:Long=0;
        Console.OUT.printf("3. Number of Aux(mu,nu) calculated at each place\nPlace  Functions  Fraction\n");
        for (i in 0..(nPlaces-1)) {
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
        Console.OUT.println ("    ROFockMethod Initialization 'up to ShellPair' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);
        
        val intPack_wlh = new WorkerLocalHandle[Integral_Pack](()=> new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        this.ylms_plh = PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs_plh();
                val ylms_plh = new Rail[Rail[Double]](shp.size); 
                finish for (i in 0..(shp.size-1)) async {
                    val sp = shp(i);
                    val tempY = new Rail[Double](sp.conA.size * sp.conB.size * (sp.maxL+1)*(sp.maxL+1));
                    intPack_wlh().genClassY(sp.bPoint, sp.aPoint, sp.zetaB, sp.zetaA, roL, tempY);
                    ylms_plh(i) = tempY;
                } 
                ylms_plh
            }
        );

        timer.stop(0);
        Console.OUT.println ("    ROFockMethod Initialization 'up to ylms' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        this.auxJMat_plh=PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs_plh(), pid = here.id;
                val auxJMat_plh = new Rail[Rail[Double]](shp.size);
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
                    auxJMat_plh(i) = new Rail[Double](size);
                }
                auxJMat_plh
            }
        );      

        this.auxKMat_plh=PlaceLocalHandle.make[Rail[Long]](
            PlaceGroup.WORLD, 
            ()=>{
                val shp = shellPairs_plh(), size=shp.size, pid=here.id;
                val auxKMat_plh = new Rail[Long](shp.size);
                finish for (i in 0..(size-1)) async {
                    val sp = shp(i), mu = sp.mu, nu = sp.nu;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu != nu) {
                        var j:Long=-1;
                        if (mu>nu) for (j=i-1; shp(j).mu!=nu || shp(j).nu!=mu; j--);
                        else for (j=i+1; shp(j).mu!=nu || shp(j).nu!=mu; j++);
                        auxKMat_plh(i)=j+size;
                        //if (shp(j).mu!=nu || shp(j).nu!=mu) Console.OUT.println ("i="+i+"j="+j);
                    } else  auxKMat_plh(i) = -1;
                }
                auxKMat_plh
            }
        );    

        val tbs = Math.max(maxRow*nOrbitals,maxam1*N)*roK; // if we use this for K too
        this.remoteBlockK_plh=PlaceLocalHandle.make[RemoteBlock](PlaceGroup.WORLD, () => new RemoteBlock(tbs/*maxRow*nOrbitals*roK*/));

        val cbs_HalfAuxInt=new Rail[Long](1); cbs_HalfAuxInt(0)=nOrbitals*roK; val halfAuxIntGrid=new Grid(cbs_HalfAuxInt, funcAtPlace);
        this.halfAuxMat=DistDenseMatrix.make(halfAuxIntGrid);

        this.auxK_wlh=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](maxam1*N*roK));
        this.dlm_wlh=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](roK));

        val cbs_nSquareMat=new Rail[Long](1); cbs_nSquareMat(0)=N; val nSquareMatGrid=new Grid(funcAtPlace, cbs_nSquareMat);
        this.distJ=DistDenseMatrix.make(nSquareMatGrid); this.distK=DistDenseMatrix.make(nSquareMatGrid);

        this.dlm_plh = PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](roK)); 
        this.e_plh = PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](2));
        this.shellPairRange_plh = shellPairRange_plh;

        this.shellAtPlace=shellAtPlace; this.offsetAtPlace=offsetAtPlace; this.funcAtPlace=funcAtPlace;
        this.roK=roK; this.roZ=jd.roZ; this.intPack_wlh=intPack_wlh; this.shellPairs_plh=shellPairs_plh;
        this.nOrbitals=nOrbitals; this.norm=bfs.getNormalizationFactors();

        timer.stop(0);
        Console.OUT.println("    ROFockMethod Initialization 'total' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        @Ifdef("__PAPI__") { papi.initialize(); papi.countFlops(); papi.countMemoryOps(); }
        @Ifdef("__DEBUG__") { printShellInfo(); }
    }

    @Native("c++", "mkl_get_max_threads()") private native static def mklGetMaxThreads():Int;
    @Native("c++", "MKL_Set_Num_Threads(#a)") private native static def mklSetNumThreads(a:Int):void;

    /** Computes the G matrix from the given density and molecular orbital matrices. */
    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}, gMatrix:DenseMatrix(N,N)) {
        Console.OUT.printf("\nROFockMethod.x10 'public def compute' %s...\n", getDateString()); 
        timer.start(TIMER_TOTAL); 

        finish ateach(place in Dist.makeUnique()) {
            val pid = here.id;
            val offsetHere = offsetAtPlace(pid);
            val nPlaces = Place.MAX_PLACES;

            ROFockMethod.setThread(1n);

            // For faster access
            val shp=shellPairs_plh(), size=shp.size, ylm=ylms_plh();
            val auxJMat=auxJMat_plh(), auxKMat=auxKMat_plh();
            val localJ=distJ.local(), localK=distK.local();
            val dlm=dlm_plh(), range=shellPairRange_plh(); 

            val bMat = new DenseMatrix(nOrbitals, roK*funcAtPlace(pid), halfAuxMat.local().d); 
            
            localJ.reset(); localK.reset();
            
            for (ron in 0n..roN) {
                // Console.OUT.println("Aux - distributed ron="+ron);
                // Aux & D
                timer.start(TIMER_AUX);
                val doK = (ron <= roNK);
                dlm.clear();
                dlm_wlh.applyLocal((d:Rail[Double])=> { d.clear(); });
                for (i in 0..(shp.size-1)) if (0<=auxKMat(i) && auxKMat(i)<size) auxKMat(i)+=size;
                
                finish DivideAndConquerLoop1D(0, shellAtPlace(pid)).execute(
                (shellIdx:Long)=> {
                    var shellPairIdx0:Long = 0; 
                    if (shellIdx>0) shellPairIdx0 = range(shellIdx-1);
                    val shellPair0 = shp(shellPairIdx0);
                    val muSize = shellPair0.mu2-shellPair0.mu+1;
                    val auxK = auxK_wlh(); 
                    auxK.clear();

                    val intPack = intPack_wlh(), dlm_partial = dlm_wlh();
                    for (var shellPairIdx:Long=shellPairIdx0; shellPairIdx<range(shellIdx); shellPairIdx++) {
                        val sp = shp(shellPairIdx);
                        val maxLron = sp.L(ron);
                        if (maxLron >= 0) {
                            val maxLm=(maxLron+1)*(maxLron+1), nuSize=sp.nu2-sp.nu+1, temp=auxJMat(shellPairIdx), srcshellPairIdx=auxKMat(shellPairIdx);
                            val off=(roK*muSize*sp.nu) as Int, asize=muSize*nuSize*roK;
                            var ind:Long=0; 
                            if (srcshellPairIdx<0 || srcshellPairIdx>=size) { // Generate Aux Ints and normalize ()
                                val y=ylm(shellPairIdx);
                                intPack.genClass(sp.bang, sp.aang, sp.bPoint, sp.aPoint, sp.zetaB, sp.zetaA, sp.conB, sp.conA, ron, maxLron, y, sp.maxL, off, auxK);
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val nrm=norm(mu)*norm(nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++) 
                                        auxK(off+ind)*=nrm;
                                }
                            } else if (doK || (srcshellPairIdx>0 && sp.nu>sp.mu)) {  // Read Aux Ints (for K or J)
                                var src:Rail[Double] = auxJMat(srcshellPairIdx); // from remote location
                                if (src.size==0) { 
                                    src=auxJMat(shellPairIdx); ind=asize;  // from local location
                                }
                                for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) 
                                    for (rolm in 0..(maxLm-1)) 
                                        auxK((nu*muSize+tmu)*roK+rolm)=src(ind++);
                            }

                            if (temp.size!=0) Rail.copy(auxK, off as Long, temp, 0, asize); // purely for J
                            else if (srcshellPairIdx>=size) { //write to remote shp J - to be read and save time
                                val temp2=auxJMat(srcshellPairIdx-size); 
                                Rail.copy(auxK, off as Long, temp2, asize, asize);
                            } 
                            if (srcshellPairIdx>=size) auxKMat(srcshellPairIdx-size)-=size;

                            // J matter
                            ind=0;
                            if (sp.mu==sp.nu) { // for diagonal block of J
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val scdmn=density(mu, nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++)
                                        dlm_partial(rolm)+=scdmn*temp(ind); 
                                }
                            } else if (temp.size!=0) { // for the rest of J
                                for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                    val scdmn=density(mu, nu);
                                    for (var rolm:Long=0; rolm<maxLm; rolm++, ind++) 
                                        dlm_partial(rolm)+=2.*scdmn*temp(ind);
                                } 
                            }
                        }
                    }
                    if (doK) {
                        val auxKMat = new DenseMatrix(roK*muSize, N, auxK);
                        DenseMatrixBLAS.compMultTrans(mos, auxKMat, bMat, [nOrbitals, auxKMat.M, N], [0, 0, 0, 0, 0, (shellPair0.mu-offsetHere)*roK], false);
                    }
                }
                );

                dlm_wlh.reduceLocal(dlm, (a:Rail[Double], b:Rail[Double])=> RailUtils.map(a, b, a, (x:Double, y:Double)=>x+y));
                Team.WORLD.allreduce[Double](dlm, 0L, dlm, 0L, dlm.size, Team.ADD);
                timer.stop(TIMER_AUX);

                finish {
                    timer.start(TIMER_K);
                    if (doK) {
                        ROFockMethod.setThread(Runtime.NTHREADS);
                        val remoteK = remoteBlockK_plh();
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
@Ifdef("__DEBUG__") {
                            val dsyrkSecs = timer.last(TIMER_DSYRK) / 1e9;
                            val dsyrkGFlops = a.N * a.N * a.M / 1e9;
                            Console.OUT.printf("Place(%d) DSYRK of %.2g GFLOPs took %.2g s ( %.2g GFLOP/s)\n", here.id, dsyrkGFlops, dsyrkSecs, (dsyrkGFlops/dsyrkSecs));
}
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

@Ifdef("__DEBUG__") {
                                val dgemmSecs = timer.last(TIMER_DGEMM) / 1e9;
                                val dgemmGFlops = 2 * a.N * thisBlock.N * a.M / 1e9;
                                Console.OUT.printf("Place(%d) DGEMM of %.2g GFLOPs from place %d took %.2g s ( %.2g GFLOP/s, %d FLOPs/word)\n", here.id, dgemmGFlops, thisBlockPlace, dgemmSecs, (dgemmGFlops/dgemmSecs), (2*a.N));
}
                            }
                        }

                        ROFockMethod.setThread(1n);
                    }
                    timer.stop(TIMER_K);
                }

                timer.start(TIMER_JMATRIX); 
                // Console.OUT.println("J - distributed"); 
                finish DivideAndConquerLoop1D(0, shp.size, 16).execute(
                (shellPairIdx:Long)=> {
                    val sp=shp(shellPairIdx);
                    val auxJ=auxJMat(shellPairIdx);
                    val maxLron=sp.L(ron);
                    if (auxJ.size > 0 && maxLron >=0n) {
                        val maxLm=(maxLron+1)*(maxLron+1);
                        var ind:Long=0; 
                        for (nu in sp.nu..sp.nu2) {
                            for (var mu:Long=sp.mu, muoff:Long=mu-offsetHere; mu<=sp.mu2; mu++, muoff++) {
                                var jContrib:Double=0.;
                                for (rolm in 0..(maxLm-1)) {
                                    jContrib += dlm(rolm) * auxJ(ind++);
                                }
                                localJ(muoff, nu) += jContrib;
                            }
                        }
                    }
                }
                );
                timer.stop(TIMER_JMATRIX);

                Team.WORLD.barrier();
            } // roN

            //Console.OUT.printf("\nG matrix\n");
            // These variables are used in the two sections below:
            val rowCount=localJ.M;
            val colCount=localJ.N;
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace((pid+1)%nPlaces);
            val colStop=offsetAtPlace((pid+mult)%nPlaces);

            // Calculate eJ and eK
            val e = e_plh();
            e.clear();
            for (var j0:Long=0, j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j0++, j++) {
                e(0) += .5*density(j, j)*localJ(j0, j);
                e(1) += .5*density(j, j)*localK(j0, j);
                for (var i0:Long=0, i:Long=offsetAtPlace(pid); i<j; i0++, i++) {
                    e(0) += density(i, j)*localJ(i0, j);
                    e(1) += density(i, j)*localK(i0, j);
                }
            }
            if (colStart < colStop) {
                for (var j:Long=colStart; j<colStop; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        e(0) +=density(ii, j)*localJ(i, j);
                        e(1) +=density(ii, j)*localK(i, j);
                    }
                }
            } else if (mult>1) {
                for (var j:Long=colStart; j<colCount; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        e(0) +=density(ii, j)*localJ(i, j);
                        e(1) +=density(ii, j)*localK(i, j);
                    }
                }
                for (var j:Long=0; j<colStop; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        e(0) +=density(ii, j)*localJ(i, j);
                        e(1) +=density(ii, j)*localK(i, j);
                    }
                }
            }
            Team.WORLD.allreduce[Double](e, 0L, e, 0L, e.size, Team.ADD);
            if (here==Place.FIRST_PLACE) Console.OUT.printf("EJ= %.10f EK= %.10f\n", e(0)/roZ, -0.5*e(1)/roZ);

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
        // gather G at place 0
        val place0GRef = new GlobalRail[Double](gMatrix.d);
        finish for (pid in 0..(nPlaces-1)) {
            val placeOffset = offsetAtPlace(pid) * gMatrix.M;
            at(Place(pid)) {
                val placeData = distJ.local().d;
                Rail.asyncCopy(placeData, 0, place0GRef, placeOffset, placeData.size);
            }
        }

        // Fix G
        for (pid in 0..(nPlaces-1)) {
            for (var j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j++) 
                for (var i:Long=offsetAtPlace(pid); i<j; i++) 
                        gMatrix(j, i) = gMatrix(i, j);
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            for (var qid:Long=pid+1; qid<pid+mult; qid++) {
                val qq=qid%nPlaces;
                for (var i:Long=offsetAtPlace(pid); i<offsetAtPlace(pid+1); i++) {
                    for (var j:Long=offsetAtPlace(qq); j<offsetAtPlace(qq+1); j++) {
                        gMatrix(j, i) = gMatrix(i, j);
                    }
                }
            }
        }
        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.2g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9); 
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
            val start = System.nanoTime();
            finish nextDim = at(Place(nextBlockPlace)) {
                val dataHere = ddm.local().d;
                Rail.asyncCopy(dataHere, 0, nextDataRef, 0, dataHere.size);
                Pair(ddm.local().M, ddm.local().N)
            };
@Ifdef("__DEBUG__") {
            val stop = System.nanoTime();
            val secs = (stop-start) / 1e9;
            val gbytes = nextDim.first * nextDim.second / 1e9;
            Console.OUT.printf("Place(%d) transfer 8 * %.2g GBs from Place(%d) took %.2g s (8 * %.2g GBytes/s)\n", here.id, gbytes, nextBlockPlace, secs, (gbytes/secs));
}
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
        @Ifdef("__MKL__") { // not working?  better use -genv OMP_NUM_THREADS 4
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
