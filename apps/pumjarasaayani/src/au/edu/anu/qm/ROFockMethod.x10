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
import x10.lang.PlaceLocalHandle;
import x10.matrix.blas.DenseMatrixBLAS;
import x10.matrix.block.Grid;
import x10.matrix.DenseMatrix;
import x10.matrix.dist.DistDenseMatrix;
import x10.regionarray.Dist;
import x10.util.Date;
import x10.util.GrowableRail;
import x10.util.Pair;
import x10.util.RailUtils;
import x10.util.Team;
import x10.util.WorkerLocalHandle;

import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ro.Integral_Pack;
import au.edu.anu.qm.ShellPair;
import au.edu.anu.util.ExecutionEnvironment;
import au.edu.anu.util.StatisticalTimer;

//import edu.utk.cs.papi.PAPI;

public class ROFockMethod(N:Long) {
    // Timer & PAPI performance counters
    public val timer=new StatisticalTimer(8);
    val TIMER_TOTAL=0;
    val TIMER_PERPLACE=1;
    val TIMER_AUX=2;
    val TIMER_JMATRIX=3;
    val TIMER_K=4;
    val TIMER_DSYRK=5;
    val TIMER_DGEMM=6;
    val TIMER_GATHER=7;

//    transient var papi:PAPI=new PAPI(); // @Ifdef("__PAPI__") // XTENLANG-3132

    val halfAuxMat:DistDenseMatrix, distJ:DistDenseMatrix, distK:DistDenseMatrix;
    val auxJ_plh:PlaceLocalHandle[Rail[Rail[Double]]], auxKIdx_plh:PlaceLocalHandle[Rail[Long]], remoteBlockK_plh:PlaceLocalHandle[RemoteBlock], ylms_plh:PlaceLocalHandle[Rail[Rail[Double]]], shellPairs_plh:PlaceLocalHandle[Rail[ShellPair]];
    val dlm_plh:PlaceLocalHandle[Rail[Double]], e_plh:PlaceLocalHandle[Rail[Double]], shellPairRange_plh:PlaceLocalHandle[Rail[Long]]; 
    val auxK_wlh:WorkerLocalHandle[Rail[Double]], intPack_wlh:WorkerLocalHandle[Integral_Pack], dlm_wlh:WorkerLocalHandle[Rail[Double]];
    val densityMos_plh:PlaceLocalHandle[Rail[DenseMatrix]];
    val nOrbitals:Long, norm:Rail[Double], roN:Int, roNK:Int, roK:Int; 
    val roZ:Double, shellAtPlace:Rail[Long], funcAtPlace:Rail[Long], offsetAtPlace:Rail[Long];

    public def this(N:Long, bfs:BasisFunctions, mol:Molecule[QMAtom], nOrbitals:Long, omega:Double, roThresh:Double) {     
        Console.OUT.printf("\nROFockMethod.x10 'public def this' %s...\n", new Date());
        property(N);
        val maxam=bfs.getShellList().getMaximumAngularMomentum();
        val maxam1=(maxam+1)*(maxam+2)/2;
        val timer=new StatisticalTimer(1), jd=JobDefaults.getInstance(), nPlaces=Place.numPlaces(),
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
        // higher angular momentum functions are assigned first.
        val nAtoms=mol.getNumberOfAtoms();
        val place2atom=new Rail[Long](nPlaces+1);
        place2atom(nPlaces)=nAtoms-1; place2atom(0)=0;
        val place2func=new Rail[Long](nPlaces+1);
        place2func(nPlaces)=mol.getAtom(nAtoms-1).getBasisFunctions().size(); place2func(0)=0;
        val npp=N/(nPlaces as Double);
        var pid:Long=nPlaces-1;
        var func:Long=0;
        for(var a:Long=nAtoms-1; a>=0 && pid>0; a--) { // centre a  
            val aFunc=mol.getAtom(a).getBasisFunctions(), naFunc=aFunc.size();          
            for(var i:Long=naFunc-1; i>=0 && pid>0; i--) { // basis functions on a
                val iaFunc=aFunc.get(i), aa=iaFunc.getTotalAngularMomentum();
                func+=(aa+1)*(aa+2)/2;
                if (N-func<=Math.ceil(pid*npp) as Long) {
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
        val roN_val=roN,roZ_val=jd.roZ,nShells_val=nShells, sizeInfo=PlaceLocalHandle.make[Rail[Double]](Place.places(), ()=> new Rail[Double](4));
        val shellPairRange_plh=PlaceLocalHandle.make[Rail[Long]](Place.places(), ()=>new Rail[Long](shellAtPlace(here.id)));
        val shellPairs_plh=PlaceLocalHandle.make[Rail[ShellPair]](
            Place.places(), 
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
        Console.OUT.printf("Fraction of N at each place: ideal=%.1f max=%.0f min=%.0f\n Imbalance cost=%.1f %%\n\n", ideal, max, min, (max/ideal-1.)*100.);
        val maxRow = max as Long; // This is used for remoteBlockK_plh allocation later

        tot=N*(1.+N)*.5; tot2=0.; max=min=at(Place(0)) sizeInfo()(0);
        var totY:Long=0, totJ:Long=0, skip:Long=0;
        Console.OUT.printf("2. Number of Aux(mu,nu) calculated at each place\nPlace  Functions  Fraction\n");
        for (i in 0..(nPlaces-1)) {
            val cost=at(Place(i)) sizeInfo()(0);
            totY+=at(Place(i)) sizeInfo()(1);
            totJ+=at(Place(i)) sizeInfo()(2);
            skip+=at(Place(i)) sizeInfo()(3);
            max=Math.max(cost,max); min=Math.min(cost,min); tot2+=cost;
            Console.OUT.printf("%5d %10.0f  %7.2f%%\n", i, cost, cost*100./tot); 
        }
        Console.OUT.printf("Fractions add up to %.2f %% (due to rounding of N/nPlaces, granularity of shellpairs and shellpair cut-off), %d shellpairs skipped\n",tot2/tot*100.0, skip);
        ideal=tot2/nPlaces;
        Console.OUT.printf("Aux/D calculations at each place: ideal=%.1f max=%.0f min=%.0f\n Imbalance cost=%.1f %%\n\n", ideal, max, min, (max/ideal-1.)*100.);

        Console.OUT.printf("Matrices size in MBs/64-bit double\nJ, K, G, density, mos\t%.3f (each)\naux4J \t%.3f\nYlm  \t%.3f\nhalfAux\t%.3f\n\n", N*N*8e-6, totJ*8e-6, totY*8e-6, nOrbitals*N*roK*8e-6);        

        timer.stop(0);
        Console.OUT.printf("    ROFockMethod Initialization 'up to ShellPair' time: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
        timer.start(0);

        val intPack_wlh = new WorkerLocalHandle[Integral_Pack](()=> new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        this.ylms_plh = PlaceLocalHandle.make[Rail[Rail[Double]]](
            Place.places(), 
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
        Console.OUT.printf("    ROFockMethod Initialization 'up to ylms' time: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
        timer.start(0);

        this.auxJ_plh=PlaceLocalHandle.make[Rail[Rail[Double]]](
            Place.places(), 
            ()=>{
                val shp = shellPairs_plh(), pid = here.id;
                val auxJ_plh = new Rail[Rail[Double]](shp.size);
                val mult = (Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                val colStart = offsetAtPlace(pid);
                val colStop = offsetAtPlace((pid+mult)%nPlaces);
                for (i in 0..(shp.size-1)) {
                    val sp = shp(i), nu = sp.nu, mu = sp.mu;
                    var size:Long = 0;
                    if ( ( (colStart<colStop) && ((colStart<=nu) && (nu<colStop)) ) ||
                         ( (colStart>=colStop) && ( (nu<colStop) || (nu>=colStart) ) ) )
                        size = sp.maxbraa*sp.maxbrab*roK;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu > nu) size=0;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu < nu) size*=2;
                    auxJ_plh(i) = new Rail[Double](size);
                }
                auxJ_plh
            }
        ); 

        timer.stop(0);
        Console.OUT.printf("    ROFockMethod Initialization 'up to auxJ' time: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
        timer.start(0);     

        this.auxKIdx_plh=PlaceLocalHandle.make[Rail[Long]](
            Place.places(), 
            ()=>{
                val shp = shellPairs_plh(), size=shp.size, pid=here.id;
                val auxKIdx_plh = new Rail[Long](shp.size);
                finish for (i in 0..(size-1)) async {
                    val sp = shp(i), mu = sp.mu, nu = sp.nu;
                    if (offsetAtPlace(pid) <= nu && nu < offsetAtPlace(pid+1) && mu != nu) {
                        var j:Long=-1;
                        if (mu>nu) for (j=i-1; shp(j).mu!=nu || shp(j).nu!=mu; j--);
                        else for (j=i+1; shp(j).mu!=nu || shp(j).nu!=mu; j++);
                        auxKIdx_plh(i)=j+size;
                        //if (shp(j).mu!=nu || shp(j).nu!=mu) Console.OUT.println ("i="+i+"j="+j);
                    } else  auxKIdx_plh(i) = -1;
                }
                auxKIdx_plh
            }
        );  

        timer.stop(0);
        Console.OUT.printf("    ROFockMethod Initialization 'up to auxK' time: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
        timer.start(0);  

        val tbs = Math.max(maxRow*nOrbitals,maxam1*N)*roK; // if we use this for K too
        this.remoteBlockK_plh=PlaceLocalHandle.make[RemoteBlock](Place.places(), () => new RemoteBlock(tbs/*maxRow*nOrbitals*roK*/));

        val cbs_HalfAuxInt=new Rail[Long](1); cbs_HalfAuxInt(0)=nOrbitals*roK; val halfAuxIntGrid=new Grid(cbs_HalfAuxInt, funcAtPlace);
        this.halfAuxMat=DistDenseMatrix.make(halfAuxIntGrid);

        this.auxK_wlh=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](maxam1*N*roK));
        this.dlm_wlh=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](roK));

        val cbs_nSquareMat=new Rail[Long](1); cbs_nSquareMat(0)=N; val nSquareMatGrid=new Grid(funcAtPlace, cbs_nSquareMat);
        this.distJ=DistDenseMatrix.make(nSquareMatGrid); this.distK=DistDenseMatrix.make(nSquareMatGrid);

        this.dlm_plh = PlaceLocalHandle.make[Rail[Double]](Place.places(), ()=> new Rail[Double](roK)); 
        this.e_plh = PlaceLocalHandle.make[Rail[Double]](Place.places(), ()=> new Rail[Double](2));
        this.shellPairRange_plh = shellPairRange_plh;

        this.shellAtPlace=shellAtPlace; this.offsetAtPlace=offsetAtPlace; this.funcAtPlace=funcAtPlace;
        this.roK=roK; this.roZ=jd.roZ; this.intPack_wlh=intPack_wlh; this.shellPairs_plh=shellPairs_plh;
        this.nOrbitals=nOrbitals; this.norm=bfs.getNormalizationFactors();

        this.densityMos_plh = PlaceLocalHandle.make[Rail[DenseMatrix]](Place.places(), ()=> new Rail[DenseMatrix](2));

        timer.stop(0);
        Console.OUT.printf("    ROFockMethod Initialization 'total' time: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
//        @Ifdef("__PAPI__") { papi.initialize(); papi.countFlops(); papi.countMemoryOps(); }
        @Ifdef("__DEBUG__") { ExecutionEnvironment.printThreadingVariables(); }
    }

    /** Computes the G matrix from the given density and molecular orbital matrices. */
    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}, gMatrix:DenseMatrix(N,N)) {
        Console.OUT.printf("\nROFockMethod.x10 'public def compute' %s...\n", new Date()); 
        timer.start(TIMER_TOTAL); 

        val tempJ = new Rail[Double](N*N);
        val tempK = new Rail[Double](N*N);
        val place0GRefJ = new GlobalRail[Double](tempJ);
        val place0GRefK = new GlobalRail[Double](tempK);

        // prepare for broadcast
        densityMos_plh()(0) = density;
        densityMos_plh()(1) = mos;

        finish ateach(place in Dist.makeUnique()) {
            timer.start(TIMER_PERPLACE);
            val nPlaces = Place.numPlaces();
            val pid = here.id;

            ExecutionEnvironment.setBlasThreads(1n);

            // For faster access
            val auxJ = auxJ_plh();
            val localJ = distJ.local(), localK = distK.local();
            val dlm = dlm_plh(); 

            val bMat = new DenseMatrix(nOrbitals, roK*funcAtPlace(pid), halfAuxMat.local().d);

            localJ.reset();
            localK.reset();
            
            if (here != Place.FIRST_PLACE) {
                densityMos_plh()(0) = new DenseMatrix(N,N);
                densityMos_plh()(1) = new DenseMatrix(N,N);
            }
            val localDensity = densityMos_plh()(0);
            val localMos = densityMos_plh()(1);
            Team.WORLD.bcast(Place.FIRST_PLACE, localDensity.d, 0, localDensity.d, 0, localDensity.d.size);
            Team.WORLD.bcast(Place.FIRST_PLACE, localMos.d, 0, localMos.d, 0, localMos.d.size);
            
            // compute contributions separately for each radial component (ron)
            for (ron in 0n..roN) {
                computeAuxBDlm(localDensity, localMos, ron, auxJ, bMat, dlm);
                computeK(ron, bMat, localK);
                computeJ(ron, auxJ, dlm, localJ);
                Team.WORLD.barrier();
            }

            timer.start(TIMER_GATHER);

            val e = e_plh();
            // gather J to tempJ, K to tempK at place 0
            finish { // only needed for timing purposes
                Rail.asyncCopy(localJ.d, 0, place0GRefJ, offsetAtPlace(pid)*N, funcAtPlace(pid)*N);
                Rail.asyncCopy(localK.d, 0, place0GRefK, offsetAtPlace(pid)*N, funcAtPlace(pid)*N);

                // while waiting for the gather to finish, calculate local contribution to energy
                // These variables are used in the two sections below:
                val rowCount = localJ.M;
                val colCount = localJ.N;
                val mult = (Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                val colStart = offsetAtPlace((pid+1)%nPlaces);
                val colStop = offsetAtPlace((pid+mult)%nPlaces);

                var eJ:Double = 0.0;
                for (var j0:Long=0, j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j0++, j++) {
                    //e(1) += .5*localDensity(j, j)*localK(j0, j);
                    for (var i0:Long=0, i:Long=offsetAtPlace(pid); i<j; i0++, i++) {
                        eJ += localDensity(i, j) * localJ(i0, j);
                        //e(1) += localDensity(i, j)*localK(i0, j);
                    }
                    eJ += 0.5 * localDensity(j, j) * localJ(j0, j);
                }
                if (colStart < colStop) {
                    for (var j:Long=colStart; j<colStop; j++) {
                        for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                            eJ += localDensity(ii, j) * localJ(i, j);
                            //e(1) +=localDensity(ii, j)*localK(i, j);
                        }
                    }
                } else if (mult > 1) {
                    for (var j:Long=colStart; j<colCount; j++) {
                        for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                            eJ += localDensity(ii, j) * localJ(i, j);
                            //e(1) +=localDensity(ii, j)*localK(i, j);
                        }
                    }
                    for (var j:Long=0; j<colStop; j++) {
                        for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                            eJ += localDensity(ii, j) * localJ(i, j);
                            //e(1) +=localDensity(ii, j)*localK(i, j);
                        }
                    }
                }
                e(0) = eJ;
            }

            Team.WORLD.allreduce[Double](e, 0L, e, 0L, e.size, Team.ADD);

            timer.stop(TIMER_GATHER);
            timer.stop(TIMER_PERPLACE);

            // Report time
            Team.WORLD.allreduce[Long](timer.total, 0L, timer.total, 0L, timer.total.size, Team.MAX);
            if (here==Place.FIRST_PLACE) {
                val tAux=(timer.total(TIMER_AUX) as Double)/1e9;
                val tJ=(timer.total(TIMER_JMATRIX) as Double)/1e9;
                val tDsyrk=(timer.total(TIMER_DSYRK) as Double)/1e9;
                val tDgemm=(timer.total(TIMER_DGEMM) as Double)/1e9;
                val tK=(timer.total(TIMER_K) as Double)/1e9;
                val tGather=(timer.total(TIMER_GATHER) as Double)/1e9;
                val tPerPlace=(timer.total(TIMER_PERPLACE) as Double)/1e9;
                Console.OUT.printf("Time (seconds) Aux= %.3f J= %.3f K-DSYRK= %.3f K-DGEMM= %.3f K-TOT= %.3f gather %.3f\n", tAux, tJ, tDsyrk, tDgemm, tK, tGather);
                Console.OUT.printf("max time per place %.3f\n", tPerPlace);
                Console.OUT.flush();
            }
        }

        val nPlaces = Place.numPlaces();
        // unpack tempK into gMatrix
        finish for (pid in 0..(nPlaces-1)) async {
            val placeRows = funcAtPlace(pid);
            val placeGOffset = offsetAtPlace(pid);
            val placeTempOffset = offsetAtPlace(pid) * N;
            for (var j:Long=0; j<N; j++)
                for (var i:Long=0; i<placeRows; i++)
                    gMatrix(i+placeGOffset, j) = -tempK(placeTempOffset + j * placeRows + i);
        }

        // Fix K to shape like J
        if (nPlaces%2==0) {
            finish for (pid in (nPlaces/2)..(nPlaces-1)) async {
                val qid=pid-nPlaces/2;
                val offset=funcAtPlace(qid)-funcAtPlace(qid)/2; // Integer rounding stuff - don't simplify
                for (var i:Long=offsetAtPlace(pid); i<offsetAtPlace(pid+1); i++)
                    for (var j:Long=offsetAtPlace(qid)+offset; j<offsetAtPlace(qid+1); j++)
                        gMatrix(i, j) = gMatrix(j, i);
            }
        }

        // add J
        finish for (pid in 0..(nPlaces-1)) async {
            val placeRows = funcAtPlace(pid);
            val placeGOffset = offsetAtPlace(pid);
            val placeTempOffset = offsetAtPlace(pid) * N;
            for (var j:Long=0; j<N; j++) 
                for (var i:Long=0; i<placeRows; i++)
                    gMatrix(i+placeGOffset, j) += tempJ(placeTempOffset + j * placeRows + i);
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
        
        val eTwo = 0.5 * traceOfSymmetricProduct(density, gMatrix) / roZ;

        val eJ = e_plh()(0)/roZ;
        val eK = (eTwo-eJ)*.5;
        timer.stop(TIMER_TOTAL);
        
        Console.OUT.printf("eTwo= %.10f EJ= %.10f EK= %.10f\n", eTwo, eJ, eK);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3f seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9); 

//        @Ifdef("__PAPI__"){ papi.printFlops(); papi.printMemoryOps();}
    }

    /** Compute trace of product of symmetric matrices a,b: trace(a &#42; b) */
    private static def traceOfSymmetricProduct(a:DenseMatrix{self.M==self.N}, b:DenseMatrix{self.M==self.N,self.M==a.M}) {
        val N = a.N;
        var trace:Double = 0.0;
        for (j in 0..(N-1)) {
            trace += a(j, j) * b(j, j);
            for (i in (j+1)..(N-1)) {
                trace += 2 * a(i, j) * b(i, j);
            }
        }
        return trace;
    }

    /** 
     * Compute auxiliary integrals, matrix B (contraction of integrals with
     * molecular orbital coefficients) and density-contracted integrals D<sub>l m</sub>.
     * @param density the density matrix Rho
     * @param mos the molecular orbital coefficient matrix C
     * @param ron the value of RO parameter n (radial component of resolution)
     * @param auxJ the auxiliary integrals in sparse-row format (output)
     * @param bMat the B matrix in dense format (output)
     * @param dlm density-contracted integrals D<sub>l m</sub> (output)
     */
    private def computeAuxBDlm(density:DenseMatrix{self.N==this.N}, mos:DenseMatrix{self.N==this.N}, ron:Int, auxJ:Rail[Rail[Double]], bMat:DenseMatrix, dlm:Rail[Double]) {
        // Console.OUT.println("Aux - distributed ron="+ron);
        timer.start(TIMER_AUX);
        //val auxStart = System.nanoTime();
        val doK = (ron <= roNK);
        dlm.clear();
        dlm_wlh.applyLocal((d:Rail[Double])=> { d.clear(); });

        val ylm = ylms_plh();
        val auxKIdx = auxKIdx_plh();

        val shp = shellPairs_plh();
        val size = shp.size;
        for (i in 0..(size-1)) if (0<=auxKIdx(i) && auxKIdx(i)<size) auxKIdx(i)+=size;
        val range = shellPairRange_plh();
        
        finish RecursiveBisection1D(0, shellAtPlace(here.id), 8).execute(
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
                    val temp = auxJ(shellPairIdx);
                    val maxLm = (maxLron+1)*(maxLron+1), nuSize = sp.nu2-sp.nu+1, srcshellPairIdx = auxKIdx(shellPairIdx);
                    val off = (roK*muSize*sp.nu) as Int, asize=muSize*nuSize*roK;
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
                        var src:Rail[Double] = auxJ(srcshellPairIdx); // from remote location
                        if (src.size==0) { 
                            src = auxJ(shellPairIdx); ind = asize;  // from local location
                        }
                        for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) 
                            for (rolm in 0..(maxLm-1)) 
                                auxK((nu*muSize+tmu)*roK+rolm) = src(ind++);
                    }

                    if (temp.size!=0) Rail.copy(auxK, off as Long, temp, 0, asize); // purely for J
                    else if (srcshellPairIdx>=size) { //write to remote shp J - to be read and save time
                        val temp2 = auxJ(srcshellPairIdx-size); 
                        Rail.copy(auxK, off as Long, temp2, asize, asize);
                    } 
                    if (srcshellPairIdx>=size) auxKIdx(srcshellPairIdx-size)-=size;

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
                DenseMatrixBLAS.compMultTrans(1.0, mos, auxKMat, 0.0, bMat, [nOrbitals, auxKMat.M, N], [0, 0, 0, 0, 0, (shellPair0.mu-offsetAtPlace(here.id))*roK]);
            }
        }
        );

        dlm_wlh.reduceLocal(dlm, (a:Rail[Double], b:Rail[Double])=> RailUtils.map(a, b, a, (x:Double, y:Double)=>x+y));
        //val auxStop = System.nanoTime();
        //Console.OUT.printf("place %d aux %d time %.3f\n", here.id, ron, (auxStop-auxStart) as Double / 1e6);
        Team.WORLD.allreduce[Double](dlm, 0L, dlm, 0L, dlm.size, Team.ADD);
        timer.stop(TIMER_AUX);
        //val allreduceStop = System.nanoTime();
        //Console.OUT.printf("place %d allreduce %d time %.3f\n", here.id, ron, (allreduceStop-auxStart) as Double / 1e6);
    }

    /** 
     * Accumulate contributions to K matrix portion local to this place:
     * K += B &#42; B<sup>T</sup>
     * @param ron the value of RO parameter n (radial component of resolution)
     * @param bMat the B matrix in dense format
     * @param localK the local portion of the K matrix
     */
    private def computeK(ron:Int, bMat:DenseMatrix, localK:DenseMatrix) {
        timer.start(TIMER_K);
        val doK = (ron <= roNK);
        if (doK) {
            finish {
                val nPlaces = Place.numPlaces();
                val pid = here.id;
                val offsetHere = offsetAtPlace(pid);
                ExecutionEnvironment.setBlasThreads(Runtime.NTHREADS);
                val remoteK = remoteBlockK_plh();
                val a = halfAuxMat.local();

                val blocks = (Math.ceil(nPlaces*.5+.5) /*- ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)*/) as Long;//
                var blk:Long = 1;
                var nextBlockPlace:Long = -1;
                finish {
                    if (blocks > 1) {
                        // overlap getting first remote block of K with DSYRK
                        nextBlockPlace = (pid+blk) % nPlaces;
                        val fullBlockSize = roK * nOrbitals * funcAtPlace(nextBlockPlace);
                        val full = (nPlaces%2==1 || blk<blocks-1 || pid<nPlaces/2);
                        val blockSize = full ? fullBlockSize : (fullBlockSize+1) / 2;
                        async remoteK.fetchNext(halfAuxMat, nextBlockPlace, blockSize);
                    }
                    timer.start(TIMER_DSYRK);
                    DenseMatrixBLAS.symRankKUpdateTrans(1.0, a, 1.0, localK, [a.N, a.M], [0, 0, 0, offsetHere], true);
                    timer.stop(TIMER_DSYRK);
@Ifdef("__DEBUG__") {
                    val dsyrkSecs = timer.last(TIMER_DSYRK) / 1e9;
                    val dsyrkGFlops = a.N * a.N * a.M / 1e9;
                    Console.OUT.printf("Place(%d) DSYRK of %.3f GFLOPs took %.3f s ( %.3f GFLOP/s)\n", here.id, dsyrkGFlops, dsyrkSecs, (dsyrkGFlops/dsyrkSecs));
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
                            val fullBlockSize = roK * nOrbitals * funcAtPlace(nextBlockPlace);
                            val full = (nPlaces%2==1 || blk<blocks-1 || pid<nPlaces/2);
                            val blockSize = full ? fullBlockSize : (fullBlockSize+1) / 2;
                            async remoteK.fetchNext(halfAuxMat, nextBlockPlace, blockSize);
                        } else {
                            nextBlockPlace = -1;
                        }
                        var dgemmGFlops:Double;
                        timer.start(TIMER_DGEMM);
                        if (nPlaces%2==1 || blk<blocks) {
                            DenseMatrixBLAS.compTransMult(1.0, a, thisBlock, 1.0, localK, [a.N, thisBlock.N, a.M], [0, 0, 0, 0, 0, offsetAtPlace(thisBlockPlace)]);
                            dgemmGFlops = 2 * a.N * thisBlock.N * a.M / 1e9;
                        } else if (pid<nPlaces/2) {
                            dgemmGFlops = 2 * (a.N/2) * thisBlock.N * a.M / 1e9;
                            DenseMatrixBLAS.compTransMult(1.0, a, thisBlock, 1.0, localK, [a.N/2, thisBlock.N, a.M], [0, a.N-a.N/2, 0, 0, a.N-a.N/2, offsetAtPlace(thisBlockPlace)]);
                        } else {
                            dgemmGFlops = 2 * (thisBlock.N-thisBlock.N/2) * thisBlock.N * a.M / 1e9;
                            DenseMatrixBLAS.compTransMult(1.0, a, thisBlock, 1.0, localK, [a.N, thisBlock.N-thisBlock.N/2, a.M], [0, 0, 0, 0, 0, offsetAtPlace(thisBlockPlace)]);
                        }
                        timer.stop(TIMER_DGEMM);

@Ifdef("__DEBUG__") {
                        val dgemmSecs = timer.last(TIMER_DGEMM) / 1e9;
                        Console.OUT.printf("Place(%d) DGEMM of %.3f GFLOPs from place %d took %.3f s ( %.3f GFLOP/s, %d FLOPs/word)\n", here.id, dgemmGFlops, thisBlockPlace, dgemmSecs, (dgemmGFlops/dgemmSecs), (2*a.N));
}
                    }
                }

                ExecutionEnvironment.setBlasThreads(1n);
            }
        }
        timer.stop(TIMER_K);
    }

    /** 
     * Accumulate contributions to J matrix portion local to this place:
     * J += D<sup>lm</sup> &#42; A<sub>mu nu, l m</sub>
     * @param ron the value of RO parameter n (radial component of resolution)
     * @param auxJ the auxiliary integrals in sparse-row format 
     * @param dlm density-contracted integrals D<sub>l m</sub> 
     * @param localJ the local portion of the J matrix
     */
    private def computeJ(ron:Int, auxJ:Rail[Rail[Double]], dlm:Rail[Double], localJ:DenseMatrix) {
        timer.start(TIMER_JMATRIX); 
        // Console.OUT.println("J - distributed");
        val shp = shellPairs_plh();
        val size = shp.size;
        val offsetHere = offsetAtPlace(here.id);
        finish RecursiveBisection1D(0, shp.size, 16).execute(
        (shellPairIdx:Long)=> {
            val sp = shp(shellPairIdx);
            val auxJSP = auxJ(shellPairIdx);
            val maxLron = sp.L(ron);
            if (auxJSP.size > 0 && maxLron >=0n) {
                val maxLm = (maxLron+1)*(maxLron+1);
                var ind:Long=0; 
                for (nu in sp.nu..sp.nu2) {
                    for (var mu:Long=sp.mu, muoff:Long=mu-offsetHere; mu<=sp.mu2; mu++, muoff++) {
                        var jContrib:Double=0.;
                        for (rolm in 0..(maxLm-1)) {
                            jContrib += dlm(rolm) * auxJSP(ind++);
                        }
                        localJ(muoff, nu) += jContrib;
                    }
                }
            }
        }
        );
        timer.stop(TIMER_JMATRIX);
    }

    private static struct RecursiveBisection1D(start:Long, end:Long, grainSize:Long) {
        public def this(start:Long, end:Long) {
            property(start, end, 1);
        }

        public def this(start:Long, end:Long, grainSize:Long) {
            property(start, end, grainSize);
        }

        public def execute(body:(idx:Long)=> void) {
            if ((end-start) > grainSize) {
                val secondHalf=RecursiveBisection1D((start+end)/2L, end, grainSize);
                async secondHalf.execute(body);
                val firstHalf=RecursiveBisection1D(start, (start+end)/2L, grainSize);
                firstHalf.execute(body);
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

        /** swap next and current */
        public def getCurrent():DenseMatrix {
            currentDim = nextDim;
            val temp = nextData;
            nextData = currentData;
            currentData = temp;
            return new DenseMatrix(currentDim.first, currentDim.second, currentData);
        }

        /** fetch next data from given place */
        public def fetchNext(ddm:DistDenseMatrix, nextBlockPlace:Long, size:Long) {
            val nextDataRef = new GlobalRail(nextData);
            val start = System.nanoTime();
            finish nextDim = at(Place(nextBlockPlace)) {
                val dataHere = ddm.local().d;
                Rail.asyncCopy(dataHere, 0, nextDataRef, 0, size);
                Pair(ddm.local().M, ddm.local().N)
            };
@Ifdef("__DEBUG__") {
            val stop = System.nanoTime();
            val secs = (stop-start) / 1e9;
            val gbytes = nextDim.first * nextDim.second / 1e9;
            Console.OUT.printf("Place(%d) transfer 8 * %.3f GBs from Place(%d) took %.3f s (8 * %.3f GBytes/s)\n", here.id, gbytes, nextBlockPlace, secs, (gbytes/secs));
}
        }

    }
}
