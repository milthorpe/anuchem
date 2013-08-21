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
    public val timer=new StatisticalTimer(4);
    val TIMER_TOTAL=0;
    val TIMER_JMATRIX=1;
    val TIMER_KMATRIX=2;
    val TIMER_GENCLASS=3;
    transient var papi:PAPI=new PAPI(); // @Ifdef("__PAPI__") // XTENLANG-3132

    // Standard conventional stuff
    private val bfs:BasisFunctions, mol:Molecule[QMAtom];
    val nOrbitals:Long, norm:Rail[Double]; 
    val shellPairs:PlaceLocalHandle[Rail[ShellPair]];

    // RO stuff 
    val auxIntMat4K:DistDenseMatrix, halfAuxMat:DistDenseMatrix, distJ:DistDenseMatrix, distK:DistDenseMatrix, tempBlock:DistDenseMatrix;
    val auxIntMat4J:PlaceLocalHandle[Rail[Rail[Double]]], ylms:PlaceLocalHandle[Rail[Rail[Double]]];
    val ttemp4J:WorkerLocalHandle[Rail[Double]], ttemp4K:WorkerLocalHandle[Rail[Double]], taux:WorkerLocalHandle[Integral_Pack], tdk:WorkerLocalHandle[Rail[Double]];
    val roZ:Double, omega:Double, roThresh:Double, funcAtPlace:Rail[Long], offsetAtPlace:Rail[Long];
    val dk:PlaceLocalHandle[Rail[Double]], e:PlaceLocalHandle[Rail[Double]];
    val roN:Int, roNK:Int, roL:Int, roK:Int; // 'var' because it can be overridden

    public def this(N:Long, bfs:BasisFunctions, mol:Molecule[QMAtom], nOrbitals:Long, omega:Double, roThresh:Double):GMatrixROmem5{self.M==N, self.N==N} {     
        super(N, N);
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def this' %s...\n", getDateString());
        val timer=new StatisticalTimer(1); timer.start(0);
        val jd=JobDefaults.getInstance(), nPlaces=Place.MAX_PLACES, nAtoms=mol.getNumberOfAtoms();
        val maxam=bfs.getShellList().getMaximumAngularMomentum(), maxam1=(maxam+1)*(maxam+2)/2;
        var nShells:Long=0;
        for (var a:Long=0; a<nAtoms; nShells+=mol.getAtom(a++).getBasisFunctions().size()) ;            

        // Set up RO roN, roL, roK for later use
        val l_n=new Rail[Int](jd.roN+3);
        val aux=new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ);
        if (omega>0.) { // long-range Ewald operator
            aux.getNL(l_n);
            roN=roNK=l_n(0);
            roL=l_n(roN+2);             
        } else { // full Coulomb operator
            roN=jd.roN;
            roL=jd.roL;
            if (jd.roNK==-1n) roNK=roN; else roNK=jd.roNK; 
        }
        val roK=(roL+1n)*(roL+1n);
        @Ifdef("__DEBUG__") {
            Console.OUT.printf("roN=%d roNK=%d roL=%d\n", roN, roNK, roL);
            Console.OUT.printf("Omega=%5.3f thresh=%e rad=%7.3f\n", omega, roThresh, jd.rad);
            Console.OUT.printf("nAtoms=%d nShells=%d N=%d\n", nAtoms, nShells, N);
            Console.OUT.printf("maxam=%d maxam1=%d\n", maxam, maxam1);
        }

        // Data distribution based on 'shell'
        // Input: nPlaces, mol
        // Output: place2atom, place2func
        val npp=N/nPlaces; var pid:Long=nPlaces-1, func:Long=0;
        val place2atom=new Rail[Long](nPlaces+1), place2func=new Rail[Long](nPlaces+1);  
        place2atom(nPlaces)=nAtoms-1; place2atom(0)=0;
        place2func(nPlaces)=mol.getAtom(nAtoms-1).getBasisFunctions().size(); place2func(0)=0;

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
            Console.ERR.println("Too many palaces!\n"); System.setExitCode(1n); 
            throw new UnsupportedOperationException("Too many places! Last pid: "+pid);
        }      
        @Ifdef("__DEBUG__") for (i in (0..(nPlaces-1))) Console.OUT.printf("Place %3d: Atom=%5d Function=%3d\n", i, place2atom(i), place2func(i));

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'Initial Assessment' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        // Generate shellpairs based on the 'shell distibution'
        // Input: mol, place2atom, place2func
        // Output: funcAtPlace, offsetAtPlace, rawShellPairs
        val funcAtPlace=new Rail[Long](nPlaces), offsetAtPlace=new Rail[Long](nPlaces+1), place2ShellPair=new Rail[Long](nPlaces+1);
        val threshold=roThresh*jd.roZ*jd.roZ*1e-3; 
        // ** Threshold must be relative to roThresh *** otherwise Z scaling will cause a problem
        // This is effectively a density threshold RO Thesis (2.26)
        val rawShellPairs=new Rail[ShellPair](nShells*nShells); // maximum limit      
        var mu:Long=0, nu:Long=0, ind:Long=0, totY:Long=0, totJ:Long=0, skip:Long=0;
        place2ShellPair(0)=0; offsetAtPlace(0)=0; offsetAtPlace(nPlaces)=N; pid=1;
        val fracJ=new Rail[Long](nPlaces);

        for(var a:Long=0; a<nAtoms; a++) { // centre a  
            val aFunc=mol.getAtom(a).getBasisFunctions();
            val naFunc=aFunc.size();          
            for (var i:Long=0; i<naFunc; i++) { // basis functions on a
                val iaFunc=aFunc.get(i);    
                val aa=iaFunc.getTotalAngularMomentum();
                val aang=iaFunc.getTotalAngularMomentum();
                val aPoint=iaFunc.origin;
                val zetaA=iaFunc.exponents;
                val conA=iaFunc.coefficients; 
                val dConA=conA.size as Int;
                val maxbraa=(aa+1)*(aa+2)/2;
                if (a==place2atom(pid) && i==place2func(pid)) {
                    place2ShellPair(pid)=ind;                    
                    funcAtPlace(pid-1)=mu-offsetAtPlace(pid-1);  
                    offsetAtPlace(pid)=mu; 
                    pid++;                      
                }          
                for(var b:Long=0; b<nAtoms; b++) { // centre b
                    val bFunc=mol.getAtom(b).getBasisFunctions();
                    val nbFunc=bFunc.size();                    
                    for(var j:Long=0; j<nbFunc; j++) { // basis functions on b
                        val jbFunc=bFunc.get(j);
                        val bb=jbFunc.getTotalAngularMomentum();                       
                        val maxbrab=(bb+1)*(bb+2)/2;     
                        val bang=jbFunc.getTotalAngularMomentum();
                        val bPoint=jbFunc.origin; 
                        val zetaB=jbFunc.exponents; 
                        val conB=jbFunc.coefficients; 
                        val dConB=conB.size as Int;
                        var contrib:Double=0.; // ss=conservative estimate
                        val R2=Math.pow(aPoint.i-bPoint.i, 2.)+Math.pow(aPoint.j-bPoint.j, 2.)+Math.pow(aPoint.k-bPoint.k, 2.);
                        for (var ii:Int=0n; ii<conA.size; ii++) {
                            for (var jj:Int=0n; jj<conB.size; jj++)  {
                                // See Szabo Ostlund 3.284-3.286
                                contrib +=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2)/Math.pow(jd.roZ, aang+bang);
                            }
                        }
                        contrib=Math.abs(contrib);
                        if (offsetAtPlace(pid-1)<=nu && nu<offsetAtPlace(pid) && mu > nu) skip++;
                        else if (contrib >=threshold) {
                            val maxL=new Rail[Int](roN+1, roL);
                            rawShellPairs(ind)=new ShellPair(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, mu, nu, maxL, contrib);     
                            ind++;
                            totY +=conA.size * conB.size * roK;
                            totJ +=maxbraa * maxbrab * roK;
                            fracJ(pid-1) += maxbraa * maxbrab; 
                        }  
                        if (b!=nAtoms-1 || j!=nbFunc-1) nu+=maxbrab; else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        place2ShellPair(pid)=ind; funcAtPlace(pid-1)=mu-offsetAtPlace(pid-1);  

        var max:Double=funcAtPlace(0), min:Double=max, tot:Double=N, ideal:Double=tot/nPlaces, tot2:Double=0;
        for (i in (1..(nPlaces-1))) {
            if (funcAtPlace(i)>max) max=funcAtPlace(i);
            if (funcAtPlace(i)<min) min=funcAtPlace(i);
        }
        val maxRow = max as Long;

        for (i in (0..(nPlaces-1))) {
            Console.OUT.printf("place %3d: offset=%5d shellpair #%5d No. of functions=%d fraction=%.2f%%\n", i, offsetAtPlace(i), place2ShellPair(i),funcAtPlace(i),funcAtPlace(i)*100./tot);
            tot2+=offsetAtPlace(i);
        }
        Console.OUT.printf("Fractions add up to %f\n",tot2/tot*100.);
        Console.OUT.printf("Load balance\n Fraction of N at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %%\n", ideal, max, min, (max/ideal-1.)*100.);

        tot=N*N*(nPlaces+1.)*.5/nPlaces;  tot2=0.;
        ideal=tot/nPlaces;
        max=min=ideal;
        for (i in (0..(nPlaces-1))) {
            pid=i;
            val mult=(Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace(pid);
            val colStop=offsetAtPlace((pid+mult)%nPlaces);
            val col=colStart<colStop?colStop-colStart:colStop+N-colStart;
            val cost=funcAtPlace(i)*col;
            Console.OUT.printf("place %3d: block size1=%d fraction=%.2f%%\n", i, cost, cost*100./tot);
            if (cost>max) max=cost;
            if (cost<min) min=cost;
            tot2+=cost;
        }
        Console.OUT.printf("Fractions add up to %f\n",tot2/tot*100.);
        Console.OUT.printf("Load balance\n Block size1 at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %%\n", ideal, max, min, (max/ideal-1.)*100.);


        tot=N*(1.+N)*.5; tot2=0.;
        ideal=tot/nPlaces;
        max=min=ideal;
        for (i in (0..(nPlaces-1))) {
            val cost=fracJ(i);
            Console.OUT.printf("place %3d: block size2=%d fraction=%.2f%%\n", i, cost, cost*100./tot);
            if (cost>max) max=cost;
            if (cost<min) min=cost;
            tot2+=cost;
        }
        Console.OUT.printf("Fractions add up to %f\n",tot2/tot*100.);
        Console.OUT.printf("Load balance\n Block size2 at each place: ideal=%.2f max=%.0f min=%.0f\n Imbalance cost=%.2f %%\n", ideal, max, min, (max/ideal-1.)*100.);


        Console.OUT.printf("Matrices larger than N-square/64-bit double/Size in MBs\naux4J \t%.3f\nYlm  \t%.3f\naux4K \t%.3f\nhalfAux\t%.3f\n", totJ*8e-6, totY*8e-6, N*N*roK*8e-6, nOrbitals*N*roK*8e-6);
        // This can be parallelised but it doesn't take significant time compare to the step below

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to ShellPair1' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        val shellPairs=PlaceLocalHandle.make[Rail[ShellPair]](
            PlaceGroup.WORLD, 
            ()=>new Rail[ShellPair](place2ShellPair(here.id+1)-place2ShellPair(here.id), 
                (i:Long)=> rawShellPairs(i+place2ShellPair(here.id))));
        // This step is very costly!

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to ShellPair2' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        val roL_val=roL;
        val taux=new WorkerLocalHandle[Integral_Pack](()=> new Integral_Pack(jd.roN, jd.roL, omega, roThresh, jd.rad, jd.roZ));
        this.ylms=PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>new Rail[Rail[Double]](place2ShellPair(here.id+1)-place2ShellPair(here.id), 
                (i:Long)=> 
                {
                    val sp=shellPairs()(i);
                    val tempY=new Rail[Double](sp.conA.size * sp.conB.size * (sp.maxL+1)*(sp.maxL+1));
                    taux().genClassY(sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, roL_val, tempY);
                    tempY
                } ));

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to ylms' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);

        this.auxIntMat4J=PlaceLocalHandle.make[Rail[Rail[Double]]](
            PlaceGroup.WORLD, 
            ()=>new Rail[Rail[Double]](place2ShellPair(here.id+1)-place2ShellPair(here.id), 
                (i:Long)=> 
                {
                    val pid=here.id, sp=shellPairs()(i);
                    val mult=(Math.ceil(nPlaces*.5+.5)-((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                    val colStart=offsetAtPlace(pid);
                    val colStop=offsetAtPlace((pid+mult)%nPlaces);
                    val nu=sp.nu; var size:Long=0L;
                    if ( ( (colStart<colStop) && ((colStart<=nu) && (nu<colStop)) ) ||
                         ( (colStart>=colStop) && ( (nu<colStop) || (nu>=colStart) ) ) )
                         size=sp.maxbraa*sp.maxbrab*roK;
                    if (offsetAtPlace(pid)<=nu && nu<offsetAtPlace(pid+1) && sp.mu > sp.nu) size=0L;
                    new Rail[Double](size)
                } ));      

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'up to auxIntMat4J' time: " + (timer.total(0) as Double) / 1e9 + " seconds");
        timer.start(0);  

        val cbs_auxInt=new Rail[Long](1); cbs_auxInt(0)=N*roK; val auxIntMat4KGrid=new Grid(funcAtPlace, cbs_auxInt);
        this.auxIntMat4K=DistDenseMatrix.make(auxIntMat4KGrid);

        val cbs_tempBlock=new Rail[Long](1); cbs_tempBlock(0)=maxRow; val tempBlockGrid=new Grid(funcAtPlace, cbs_tempBlock);
        this.tempBlock=DistDenseMatrix.make(tempBlockGrid);

        val cbs_HalfAuxInt=new Rail[Long](1); cbs_HalfAuxInt(0)=nOrbitals*roK; val halfAuxIntGrid=new Grid(funcAtPlace, cbs_HalfAuxInt);
        this.halfAuxMat=DistDenseMatrix.make(halfAuxIntGrid);

        this.ttemp4J=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](maxam1*maxam1*roK));
        this.ttemp4K=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](maxam1*maxam1*roK));
        this.tdk=new WorkerLocalHandle[Rail[Double]](()=> new Rail[Double](roK));

        val cbs_nSquareMat=new Rail[Long](1); cbs_nSquareMat(0)=N; val nSquareMatGrid=new Grid(funcAtPlace, cbs_nSquareMat);
        this.distJ=DistDenseMatrix.make(nSquareMatGrid); this.distK=DistDenseMatrix.make(nSquareMatGrid);

        this.dk=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](roK)); 
        this.e=PlaceLocalHandle.make[Rail[Double]](PlaceGroup.WORLD, ()=> new Rail[Double](2));
        this.norm=bfs.getNormalizationFactors(); // Vector quantity

        this.roK=roK; this.roZ=jd.roZ;
        this.taux=taux; this.shellPairs=shellPairs;
        this.offsetAtPlace=offsetAtPlace; this.funcAtPlace=funcAtPlace;
        this.bfs=bfs; this.mol=mol; this.nOrbitals=nOrbitals; this.omega=omega; this.roThresh=roThresh;

        timer.stop(0);
        Console.OUT.println ("    GMatrixROmem5 Initialization 'total' time: " + (timer.total(0) as Double) / 1e9 + " seconds");

        @Ifdef("__MKL__") {
            Console.OUT.print("mklGetMaxThreads() was " + mklGetMaxThreads() + " and is now set to"); mklSenumThreads(Runtime.NTHREADS);
            Console.OUT.println(" " + mklGetMaxThreads() + " thread(s).");
            // not working on multiple nodes?  better use -env OMP_NUM_THREADS 4
        }  
        @Ifdef("__PAPI__") { papi.initialize(); papi.countFlops(); papi.countMemoryOps();}
    }

    @Native("c++", "mkl_get_max_threads()") private native static def mklGetMaxThreads():Int;
    @Native("c++", "MKL_Set_Num_Threads(#a)") private native static def mklSenumThreads(a:Int):void;

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        Console.OUT.printf("\nGMatrixROmem5.x10 'public def compute' %s...\n", getDateString()); 

        @Ifdef("__DEBUG__") finish ateach(place in Dist.makeUnique()) checkTopology();            
        // Console.OUT.flush();

        timer.start(TIMER_TOTAL); 

        val timer=this.timer; 
        val shellPairs=this.shellPairs; val ylms=this.ylms, distJ=this.distJ, distK=this.distK; 
        val N=this.N, nOrbitals=this.nOrbitals, roN=this.roN, roNK=this.roNK, roL=this.roL, roK=this.roK, roZ=this.roZ, norm=this.norm, e=this.e;
        val funcAtPlace=this.funcAtPlace, offsetAtPlace=this.offsetAtPlace, taux=this.taux, dk=this.dk;
        val auxIntMat4K=this.auxIntMat4K, halfAuxMat=this.halfAuxMat, tempBlock=this.tempBlock;
        val auxIntMat4J=this.auxIntMat4J, ttemp4J=this.ttemp4J, ttemp4K=this.ttemp4K, tdk=this.tdk;
      
        finish ateach(place in Dist.makeUnique()) {
            val pid=here.id;
            val nPlaces=Place.MAX_PLACES;

            // For faster access 
            val shp=shellPairs(), ylmp=ylms();
            val localAuxJ=auxIntMat4J(), localMatK=auxIntMat4K.local(), tBlock=tempBlock.local(); 
            val localJ=distJ.local(), localK=distK.local();
            val dkp=dk(), ep=e(); 
            
            ep.clear(); localJ.reset(); localK.reset();

            for (ron in 0n..roN) {
                // Console.OUT.println("Aux - distributed ron="+ron);
                timer.start(TIMER_GENCLASS); 
                dkp.clear(); tdk.applyLocal((d:Rail[Double])=> { d.clear(); });

                // Aux & D
                finish DivideAndConquerLoop1D(0, shp.size).execute(
                (spInd:Long)=> {
                    val sp=shp(spInd);
                    val maxLron=sp.L(ron);
                    if (maxLron >=0) {
                        var auxJ:Rail[Double]=localAuxJ(spInd);
                        if (auxJ.size==0L) auxJ=ttemp4J();
                        val aux=taux(), auxK=ttemp4K(), myThreaddk=tdk();
                        val maxLm=(maxLron+1)*(maxLron+1), y=ylmp(spInd);
                        val musize=sp.mu2-sp.mu+1, nusize=sp.nu2-sp.nu+1;
                        aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, ron, maxLron, y, sp.maxL, auxJ);
                        var ind:Long=0;
                        if (localAuxJ(spInd).size==0L) {
                            for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) {
                                for (var nu:Long=sp.nu, tnu:Long=0; nu<=sp.nu2; nu++, tnu++) {
                                    val nrm=norm(mu)*norm(nu);
                                    for (rolm in 0..(maxLm-1)) {
                                        val normAux=(auxJ(ind++) *=nrm);
                                        auxK((tnu*roK+rolm)*musize+tmu)=normAux;
                                    }
                                }
                            }
                        } else if (sp.mu!=sp.nu) {
                            for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) for (var nu:Long=sp.nu, tnu:Long=0; nu<=sp.nu2; nu++, tnu++) {
                                val scdmn=density(mu, nu);
                                val nrm=norm(mu) * norm(nu); 
                                for (rolm in 0..(maxLm-1)) {
                                    val normAux=(auxJ(ind++) *=nrm);
                                    myThreaddk(rolm) +=2.*scdmn*normAux; 
                                    auxK((tnu*roK+rolm)*musize+tmu)=normAux;
                                }
                            }
                        } else {
                            for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) {
                                for (var nu:Long=sp.nu, tnu:Long=0; nu<=sp.nu2; nu++, tnu++) {
                                    val scdmn=density(mu, nu);
                                    val nrm=norm(mu) * norm(nu); 
                                    for (rolm in 0..(maxLm-1)) {
                                        val normAux=(auxJ(ind++) *=nrm);
                                        myThreaddk(rolm) +=scdmn*normAux; 
                                        auxK((tnu*roK+rolm)*musize+tmu)=normAux;
                                    }
                                }
                            }
                        }
                        ind=0;
                        val rows=sp.mu2-sp.mu+1;
                        for (var nu:Long=sp.nu; nu<=sp.nu2; nu++) {
                            DenseMatrix.copySubset(auxK, ind, localMatK, sp.mu-offsetAtPlace(pid), nu*roK, rows, maxLm);
                            ind +=rows*maxLm;
                        }
                        if (offsetAtPlace(pid) <=sp.nu && sp.nu < offsetAtPlace(pid+1) && sp.mu < sp.nu) {
                            ind=0;
                            for (var mu:Long=sp.mu, tmu:Long=0; mu<=sp.mu2; mu++, tmu++) {
                                for (var nu:Long=sp.nu, tnu:Long=0; nu<=sp.nu2; nu++, tnu++) {
                                    for (rolm in 0..(maxLm-1)) {
                                        auxK((tmu*roK+rolm)*nusize+tnu)=auxJ(ind++);
                                    }
                                }
                            }
                            val rows2=sp.nu2-sp.nu+1;
                            ind=0;
                            for (var mu:Long=sp.mu; mu<=sp.mu2; mu++) {
                                DenseMatrix.copySubset(auxK, ind, localMatK, sp.nu-offsetAtPlace(pid), mu*roK, rows2, maxLm);
                                ind +=rows2*maxLm;
                            }
                        }
                    }
                }
                );

                tdk.reduceLocal(dkp, (a:Rail[Double], b:Rail[Double])=> RailUtils.map(a, b, a, (x:Double, y:Double)=>x+y));
                timer.stop(TIMER_GENCLASS);

                finish {
                    async Team.WORLD.allreduce[Double](dkp, 0L, dkp, 0L, dkp.size, Team.ADD);
                    timer.start(TIMER_KMATRIX);
                    if (ron <=roNK) {
                        val A=new DenseMatrix(funcAtPlace(pid)*roK, N, localMatK.d);
                        val B=new DenseMatrix(funcAtPlace(pid)*roK, nOrbitals, halfAuxMat.local().d); 
                        DenseMatrixBLAS.compMultTrans(A, mos, B, [funcAtPlace(pid)*roK, nOrbitals, N], false);
                        //cannot do B.multTrans(A, mos, false); -  mos is [N, N] rather than [nObital, N]
                    }
                    timer.stop(TIMER_KMATRIX);
                }

                if (ron <=roNK) { // This produces K/2 & can be improved further by better scheduling and buffering
                    timer.start(TIMER_KMATRIX);
                    // Console.OUT.println("K - distributed");
                    val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
                    val a=halfAuxMat.local();
                    for (var blk:Long=0; blk<mult; blk++) {
                        val qid=(pid+blk)%nPlaces;
                        val b=at(Place(qid)) {halfAuxMat.local()}; // This might be a waste of memory?
                        val noff=offsetAtPlace(qid);
                        val c=new DenseMatrix(a.M, b.M, tBlock.d); 
                        c.multTrans(a, b, false);
                        for (var j:Long=0; j<c.N; j++) for (var i:Long=0; i<c.M; i++) {
                            localK(i, noff+j) +=c(i, j);
                        }
                    }
                    timer.stop(TIMER_KMATRIX);
                }

                timer.start(TIMER_JMATRIX); 
                // Console.OUT.println("J - distributed"); 
                finish DivideAndConquerLoop1D(0, shp.size).execute(
                (spInd:Long)=> {
                    val sp=shp(spInd);
                    val auxJ=localAuxJ(spInd);
                    val maxLron=sp.L(ron);
                    if (auxJ.size > 0 && maxLron >=0n) {
                        val maxLm=(maxLron+1)*(maxLron+1);
                        var ind:Long=0; 
                        for (var mu:Long=sp.mu, muoff:Long=mu-offsetAtPlace(pid); mu<=sp.mu2; mu++, muoff++) {
                            for (nu in sp.nu..sp.nu2) {
                                var jContrib:Double=0.;
                                for (rolm in 0..(maxLm-1)) {
                                    jContrib +=dkp(rolm) * auxJ(ind++);
                                }
                                localJ(muoff, nu) +=jContrib;
                            }
                        }
                    }
                }
                );
                timer.stop(TIMER_JMATRIX);
                Team.WORLD.barrier();
            }

            //Console.OUT.printf("\nG matrix\n");

            // Fix J
            finish DivideAndConquerLoop1D(0, shp.size).execute(
            (spInd:Long)=> {
                val sp=shp(spInd);
                if (localAuxJ(spInd).size > 0 && offsetAtPlace(pid) <=sp.nu && sp.nu < offsetAtPlace(pid+1) && sp.mu !=sp.nu) {
                    for (var mu:Long=sp.mu, muoff:Long=mu-offsetAtPlace(pid); mu<=sp.mu2; mu++, muoff++) {
                        for (var nu:Long=sp.nu, nuoff:Long=nu-offsetAtPlace(pid); nu<=sp.nu2; nu++, nuoff++) {
                            localJ(nuoff, mu)=localJ(muoff, nu);
                        }
                    } 
                }
            }
            );

            // These variables are used in the two sections below:
            val rowCount=localJ.M;
            val colCount=localJ.N;
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            val colStart=offsetAtPlace(pid);
            val colStop=offsetAtPlace((pid+mult)%nPlaces);

            // Calculate eJ and eK
            for (var j:Long=offsetAtPlace(pid); j<offsetAtPlace(pid+1); j++) {
                for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                    ep(0) -=.5*density(ii, j)*localJ(i, j);
                    ep(1) -=.5*density(ii, j)*localK(i, j);
                }
            }

            if (colStart < colStop) {
                for (var j:Long=colStart; j<colStop; j++) {
                    for (var i:Long=0, ii:Long=offsetAtPlace(pid); i<rowCount; i++, ii++) {
                        ep(0) +=density(ii, j)*localJ(i, j);
                        ep(1) +=density(ii, j)*localK(i, j);
                    }
                }
            } else {
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
            if (colStart < colStop) {
                for (var j:Long=colStart; j<colStop; j++) {
                    for (var i:Long=0; i<rowCount; i++) {
                        localJ(i, j) -=localK(i, j);
                    }
                }
            } else {
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
                val tINT=(timer.total(TIMER_GENCLASS) as Double)/1e9;
                val tJ=(timer.total(TIMER_JMATRIX) as Double)/1e9;
                val tK=(timer.total(TIMER_KMATRIX) as Double)/1e9;
                Console.OUT.printf("Time INT= %.2f s J= %.2f s K= %.2f s\n", tINT, tJ, tK);
                Console.OUT.flush();
            }
        }

        val nPlaces=Place.MAX_PLACES;
        // Copy G to place 0 // can improve further by copying only contributing blocks
        for (pid in 0..(nPlaces-1)) {
            val mat=at(Place(pid)) { distJ.local() };
            DenseMatrix.copySubset(mat, 0, 0, this, offsetAtPlace(pid), 0, mat.M, mat.N);
        }

        // Fix G
        for (pid in 0..(nPlaces-1)) {
            val mult=(Math.ceil(nPlaces*.5+.5) - ((nPlaces%2L==0L && pid<nPlaces/2)?1:0)) as Long;
            for (var qid:Long=mult+pid; qid<nPlaces+pid; qid++) {
                val qq=qid%nPlaces;
                for (var i:Long=offsetAtPlace(pid); i<offsetAtPlace(pid+1); i++) {
                    for (var j:Long=offsetAtPlace(qq); j<offsetAtPlace(qq+1); j++) {
                        this(i, j)=this(j, i);
                    }
                }
            }
        }
        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9); 
        @Ifdef("__PAPI__"){ papi.printFlops(); papi.printMemoryOps();}
    }

    private static struct DivideAndConquerLoop1D(start:Long, end:Long) {
        // TODO allow grain size > 1
        public def this(start:Long, end:Long) {
            property(start, end);
        }

        public def execute(body:(idx:Long)=> void) {
            if ((end-start) > 1L) {
                val firstHalf=DivideAndConquerLoop1D(start, (start+end)/2L);
                val secondHalf=DivideAndConquerLoop1D((start+end)/2L, end);
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
            val dateReader=Runtime.execForRead("date"); 
            result=dateReader.readLine();
        } catch (e:IOException) {
            // could not read date! use current time in milliseconds
            result=System.currentTimeMillis() + "ms";
        }
        return result;
    }

    private def checkTopology() {
        val hostnameReader=Runtime.execForRead("uname -n"); 
        val hostname=hostnameReader.readLine();
        
        val npReader=Runtime.execForRead("echo $X10_NPLACES");
        val np=npReader.readLine();

        val ntReader=Runtime.execForRead("echo $X10_NTHREADS");
        val nt=ntReader.readLine();

        val gtReader=Runtime.execForRead("echo $GOTO_NUM_THREADS");
        val gt=gtReader.readLine();

        val ompReader=Runtime.execForRead("echo $OMP_NUM_THREADS");
        val omp=ompReader.readLine();

        Console.OUT.println(here + ", Runtime.NTHREADS=" + Runtime.NTHREADS + ", uname -n=" + hostname + ", X10_NPLACES="+np+", X10_NTHREADS="+nt+", GOTO_NUM_THREADS="+gt+ ", OMP_NUM_THREADS="+omp ); // if print out separately, it goes randomly.
        Console.OUT.flush();
    }
}
