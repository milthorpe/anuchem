/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;

import x10.util.ArrayList;
import x10.util.ArrayUtils;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;

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

public class GMatrixROmem3 extends DenseMatrix{self.M==self.N} {
    // Timer & PAPI performance counters
    public val timer = new StatisticalTimer(4);
    static TIMER_TOTAL = 0; static TIMER_JMATRIX = 1;
    static TIMER_KMATRIX = 2; static TIMER_GENCLASS = 3;
    // @Ifdef("__PAPI__") // XTENLANG-3132
    transient var papi:PAPI = new PAPI(); 

    // RO stuff 
    transient val aux:Integral_Pack;
    var roN:Int; var roNK:Int; var roL:Int; // 'var' because it can be overridden
    val roK:Int; val roZ:Double; val temp:Rail[Double]; 
    val auxIntMat:DenseMatrix; 
    val halfAuxMat:DenseMatrix; 
    val ylms:Rail[Ylm];

    // Standard conventional stuff
    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
    val nOrbital:Int; val numSigShellPairs:Int;
    val norm:Rail[Double]; val dk:Rail[Double];
    val emptyRailD = new Rail[Double](0), emptyRailI = new Rail[Int](0); 
    val shellPairs:Rail[ShellPair]; 
    val jMatrix:DenseMatrix{self.M==self.N,self.N==this.N};
    val kMatrix:DenseMatrix{self.M==self.N,self.N==this.N};

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int, omega:Double,roThresh:Double):GMatrixROmem3{self.M==N,self.N==N} {     
        super(N, N); 
        val result = Runtime.execForRead("date"); Console.OUT.printf("\nGMatrixROmem.x10 'public def this' %s...\n",result.readLine()); 
        this.bfs = bfs; this.mol = molecule; this.nOrbital = nOrbital;     

        @Ifdef("__PAPI__") { 
            papi.initialize();
            papi.countFlops();
            papi.countMemoryOps(); 
        }

        val jd = JobDefaults.getInstance();
        val l_n = new Rail[Int](jd.roN+3);
        aux = new Integral_Pack(jd.roN,jd.roL,omega,roThresh,jd.rad);          
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
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roK); // will be passed to C++ code           
        dk = new Rail[Double](roK); // eqn 15b in RO#7
        var maxint:Rail[Double] = new Rail[Double](roK*(roN+1)); 
     
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N, N);
        kMatrix = new DenseMatrix(N, N);
        @Ifdef("__DEBUG__") { Console.OUT.printf("GMatrixROmem3.x10 Memory O(N^2) allocated sucessfully...\n"); }

        auxIntMat = new DenseMatrix(N*roK,N);
        halfAuxMat = new DenseMatrix(N*roK,nOrbital);
        @Ifdef("__DEBUG__") { Console.OUT.printf("GMatrixROmem3.x10 Memory O(N^2 K) allocated sucessfully...\n"); }        

        // Shell/Shellpair business 
        Console.OUT.printf("GMatrixROmem.x10 Screening shellpairs...\n");        
        val threshold = roThresh*jd.roZ*jd.roZ*1e-3; // ** must be relative to roThresh *** otherwise Z scaling will cause a problem
        // This is effectively a density threshold RO Thesis (2.26)
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }
        var rawShellPairs:Rail[ShellPair] = new Array[ShellPair](nShell*(nShell+1)/2); 
        val zeroPoint = Point3d(0.0, 0.0, 0.0);
       
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRailD, emptyRailD, emptyRailD, emptyRailD, 0, 0, 0, 0, emptyRailI, threshold);
        var mu:Int=0,nu:Int=0,ind:Int=0; 
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
                        //Console.OUT.printf("a=%d i=%d b=%d j=%d [naFunc=%d nbFunc=%d]\n", a,i,b,j,naFunc,nbFunc);
                        var aaFunc:ContractedGaussian=iaFunc,bbFunc:ContractedGaussian=jbFunc;
                        val aa=iaFunc.getTotalAngularMomentum(); val bb=jbFunc.getTotalAngularMomentum();                       
                        val maxbraa = (aa+1)*(aa+2)/2; val maxbrab = (bb+1)*(bb+2)/2;                                          
                        if (aa>bb || (aa==bb && mu>=nu)) { // Careful this is a tricky condition                            
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
                            }               
                        }
 
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab;
                        else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        numSigShellPairs = ind;
        Console.OUT.printf("Found %d significant shellpairs.\n",numSigShellPairs);
        shellPairs = new Rail[ShellPair](numSigShellPairs);    
        ylms = new Rail[Ylm](numSigShellPairs);     
        for (i in 0..(numSigShellPairs-1)) {
            val sh=rawShellPairs(i); 
            shellPairs(i)=sh;            
            val tempY=new Rail[Double](sh.dconA*sh.dconB*(roL+1)*(roL+1));
            aux.genClassY(sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.dconA, sh.dconB, roL, tempY);
            ylms(i) = new Ylm(tempY,roL);
        }
        rawShellPairs = null; // Deallocate this variable     
    }

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        val result = Runtime.execForRead("date"); Console.OUT.printf("\nGMatrixROmem.x10 'public def compute' %s...\n",result.readLine()); 
        timer.start(TIMER_TOTAL); var tINT:Double=0.,tJ:Double=0.,tK:Double=0.;
        jMatrix.reset(); kMatrix.reset(); var ron:Int=0;         
        /*finish*/ /*for (n in 0..roN)*/for (; ron<=roN; ron++) /*async*/ {
            @Ifdef("__DEBUG__") {Console.OUT.printf("ron=%d...\n",ron); }
            dk.clear();     

            timer.start(TIMER_GENCLASS);   
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd); val maxLron=sp.maxL(ron);                     
                if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1); var ind:Int=0;                                                              
                    aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron,ylms(spInd).y, ylms(spInd).maxL);
                    if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++)/*for ([tmu,tnu] in (sp.mu..sp.mu2)*(sp.nu..sp.nu2))*/ {
                        val scdmn=density(tmu,tnu)*2.;
                        val nrm=norm(tmu)*norm(tnu);
                        for (var rolm:Int=0; rolm<maxLm; rolm++) {
                            val normAux = nrm*temp(ind++); 
                            /*atomic*/ dk(rolm) += scdmn*normAux; 
                            auxIntMat(tmu+rolm*N,tnu) = auxIntMat(tnu+rolm*N,tmu) = normAux;                            
                        } 
                    }
                    else for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) /*for ([tmu,tnu] in (sp.mu..sp.mu2)*(sp.nu..sp.nu2))*/ {
                        val scdmn=density(tmu,tnu);
                        val nrm=norm(tmu)*norm(tnu);
                        for (var rolm:Int=0; rolm<maxLm; rolm++) {
                            val normAux = nrm*temp(ind++);
                            /*atomic*/ dk(rolm) += scdmn*normAux; 
                            auxIntMat(tmu+rolm*N,tnu) = normAux;                            
                        } 
                    }

                }
            }
            timer.stop(TIMER_GENCLASS); tINT+=(timer.last(TIMER_GENCLASS) as Double)/1e9;

            timer.start(TIMER_JMATRIX);             
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (sp.maxL(ron)>=0) { 
                    val maxLm=(maxLron+1)*(maxLron+1); 
                    for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) /*for ([tmu,tnu] in (sp.mu..sp.mu2)*(sp.nu..sp.nu2))*/ {
                        var jContrib:Double=0.;
                        for (var rolm:Int=0; rolm<maxLm; rolm++) jContrib+=dk(rolm)*auxIntMat(tmu+rolm*N,tnu);
                        /*atomic*/ jMatrix(tmu,tnu) += jContrib;
                    } 
                }
            }            
            timer.stop(TIMER_JMATRIX); tJ+=(timer.last(TIMER_JMATRIX) as Double)/1e9;
        
            if (ron<=roNK) { //This produces K/2
                timer.start(TIMER_KMATRIX);  
                /*atomic*/ DenseMatrixBLAS.compMultTrans(auxIntMat, mos, halfAuxMat, [N*roK, nOrbital, N], false);
                val halfAuxMat2 = new DenseMatrix(N, roK*nOrbital, halfAuxMat.d);
                /*atomic*/ DenseMatrixBLAS.compMultTrans(halfAuxMat2, halfAuxMat2, kMatrix, [N, N, roK*nOrbital], true);
                timer.stop(TIMER_KMATRIX); tK+=(timer.last(TIMER_KMATRIX) as Double)/1e9;
            }
        }

        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }

        // Fix upper half of J (see the definition of shellPairs) ==> sp.mu>sp.nu
        // Fix the whole matrix ==> if (sp.mu!=sp.nu)
        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<=sp.mu2; tmu++) for (var tnu:Int=sp.nu; tnu<=sp.nu2; tnu++) /*for ([tmu,tnu] in (sp.mu..sp.mu2)*(sp.nu..sp.nu2))*/
                jMatrix(tnu,tmu) = jMatrix(tmu,tnu);
        }
        // Use the upper half of J and K to form G
        for (var tmu:Int=0; tmu<N; tmu++) for (var tnu:Int=tmu; tnu<N; tnu++) 
            this(tnu,tmu)=this(tmu,tnu)=jMatrix(tmu,tnu)-kMatrix(tmu,tnu);

        @Ifdef("__DEBUG__") {
            var eJ:Double=0.,eK:Double=0.;
            for (var tmu:Int=0; tmu<N; tmu++) for (var tnu:Int=tmu+1; tnu<N; tnu++) {
                eJ+=density(tmu,tnu)*jMatrix(tmu,tnu);
                eK-=.5*density(tmu,tnu)*kMatrix(tmu,tnu);
            }
            for (var tmu:Int=0; tmu<N; tmu++) {
                eJ+=.5*density(tmu,tmu)*jMatrix(tmu,tmu);
                eK-=.25*density(tmu,tmu)*kMatrix(tmu,tmu);
            }
            Console.OUT.printf("  EJ = %.10f a.u.\n  EK = %.10f a.u.\n", eJ/roZ, eK/roZ);
            Console.OUT.printf("  Time INT = %.2f s J = %.2f s K = %.2f s\n", tINT, tJ, tK);
        }

        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);
    }
}
