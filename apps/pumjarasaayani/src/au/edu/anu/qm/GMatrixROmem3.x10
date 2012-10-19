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
    // Timer
    public val timer = new StatisticalTimer(7);
    static TIMER_TOTAL = 0; static TIMER_JMATRIX = 1;
    static TIMER_KMATRIX = 2; static TIMER_GENCLASS = 3;

    static TIMER_INTGEN = 4;
    static TIMER_INTSCREEN = 5;
    static TIMER_INTCOPY = 6;
    // PAPI performance counters
    // @Ifdef("__PAPI__") // XTENLANG-3132
    transient var papi:PAPI = new PAPI(); 

    // Integral_Pack 
    transient val aux:Integral_Pack;
    var roN:Int; var roNK:Int; var roL:Int; // 'var' because it can be overridden
    val roK:Int; val temp:Rail[Double]; 

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
    val nOrbital:Int; val numSigShellPairs:Int;
    val norm:Rail[Double]; val dk:Rail[Double];
    val emptyRailD = new Rail[Double](0), emptyRailI = new Rail[Int](0); 
    val shellPairs:Rail[ShellPair]; 
    val ylms:Rail[Ylm];

    val jMatrix:DenseMatrix{self.M==self.N,self.N==this.N};
    val kMatrix:DenseMatrix{self.M==self.N,self.N==this.N};
    val scratch:DenseMatrix{self.M==self.N,self.N==this.N};

    // Big arrays stored in Memory     
    val auxIntMat:DenseMatrix; 
    val halfAuxMat:DenseMatrix; 

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
        aux = new Integral_Pack(jd.roN,jd.roL,omega,roThresh,jd.rad*jd.roZ);        
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
        roK = (roL+1)*(roL+1);        

        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val mdc=bfs.getShellList().getMaximumDegreeOfContraction();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roK); // will be passed to C++ code           
        dk = new Rail[Double](roK); // eqn 15b in RO#7
        var maxint:Rail[Double] = new Rail[Double](roK*(roN+1)); 
     
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N, N);
        kMatrix = new DenseMatrix(N, N);
        scratch = new DenseMatrix(N, N);
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
        var rawShellPairs:Array[ShellPair] = new Array[ShellPair](nShell*(nShell+1)/2); 
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
                        if (aa>bb || (aa==bb && mu>=nu)) {
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
        jMatrix.reset(); kMatrix.reset();
        val jd = JobDefaults.getInstance();     
        
        var fac:Double; var ind:Int; var jContrib:Double;
        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) fac=2.; else fac=1.;
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++)
                scratch(tmu,tnu)=density(tmu,tnu)*fac*norm(tmu)*norm(tnu);            
        }

        for (var ron:Int=0; ron<=roN; ron++) {
            @Ifdef("__DEBUG__") {Console.OUT.printf("ron=%d...\n",ron); }
            dk.clear();       
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);                
                if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1); ind=0;                                         
                    aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron,ylms(spInd).y, ylms(spInd).maxL);
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                        val scdmn=scratch(tmu,tnu);
                        val muOffset=tmu*roK;
                        val nrm=norm(tmu)*norm(tnu);
                        for (var rolm:Int=0; rolm<maxLm; rolm++) {
                            dk(rolm) += scdmn*temp(ind); 
                            auxIntMat(muOffset+rolm,tnu) = nrm*temp(ind);
                            ind++;  
                        } 
                        if (sp.mu!=sp.nu) {
                            val nuOffset=tnu*roK;
                            for (var rolm:Int=0; rolm<maxLm; rolm++) auxIntMat(nuOffset+rolm,tmu)=auxIntMat(muOffset+rolm,tnu);
                        }
                    }
                }
            }
             
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (sp.maxL(ron)>=0) { 
                    val maxLm=(maxLron+1)*(maxLron+1); ind=0;
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                        val nrm=norm(tmu)*norm(tnu);
                        jContrib = 0.;
                        val muOffset=tmu*roK;
                        for (var rolm:Int=0; rolm<maxLm; rolm++)
                            jContrib+=dk(rolm)*nrm*auxIntMat(muOffset+rolm,tnu);
                        jMatrix(tmu,tnu) += jContrib;
                    } 
                }
            }            
                    
            if (ron<=roNK) {        
                for (var i:Int=0; i<N*roK; i++) for (var j:Int=0; j<nOrbital; j++) {
                    var contrib:Double=0.;
                    for (var k:Int=0; k<N; k++) contrib+= mos(j,k)*auxIntMat(i,k);
                    halfAuxMat(i,j)=contrib;
                }
                                                
                for (var tmu:Int=0; tmu<N; tmu++) for (var tnu:Int=0; tnu<N; tnu++) {
                    var kContrib:Double=0.;
                    val muOffset=tmu*roK; val nuOffset=tnu*roK;
                    for (var orb:Int=0; orb<nOrbital; orb++) for (var rolm:Int=0; rolm<roK; rolm++)
                        kContrib+=halfAuxMat(muOffset+rolm,orb)*halfAuxMat(nuOffset+rolm,orb);
                    kMatrix(tmu,tnu) += kContrib;           
                }
            }
        }

        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) 
                jMatrix(tnu,tmu) = jMatrix(tmu,tnu);
        }      
        
        @Ifdef("__PAPI__"){
            papi.printFlops();
            papi.printMemoryOps();
        }

        @Ifdef("__DEBUG__") {
            val eJ = scratch.mult(density, jMatrix).trace();
            Console.OUT.printf("  EJ = %.10f a.u.\n", eJ/jd.roZ);
            val eK = scratch.mult(density, kMatrix).trace();
            Console.OUT.printf("  EK = %.10f a.u.\n", eK/jd.roZ);
        }
        // Form G matrix
        jMatrix.d.map(this.d, kMatrix.d, (j:Double,k:Double)=>(2.0*j-k)); // eqn 14
        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);
    }
}
