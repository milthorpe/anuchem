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

import x10.util.ArrayList;
import x10.util.ArrayUtils;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;

import x10x.vector.Vector;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ShellPair; 
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;

/**
 * G matrix in HF calculation -- RO 
 */

public class GMatrixROmem extends DenseMatrix{self.M==self.N} {

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    val roN:Int;
    val roNK:Int;
    val roL:Int;
    val roZ:Double;
    val nOrbital:Int;
    
    val norm:Rail[Double];
    val temp:Rail[Double];
    val dk:Rail[Double];
    val muk:DenseMatrix{self.M==this.N};

    val shellPairs:Rail[ShellPair];
    val numSigShellPairs:Int;

    transient val aux:Integral_Pack;
    val jMatrix:DenseMatrix{self.M==self.N,self.N==this.N};
    val kMatrix:DenseMatrix{self.M==self.N,self.N==this.N};

    // used for calculating eJ, eK
    val scratch:DenseMatrix{self.M==self.N,self.N==this.N};

    var counter:Int=0;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int):GMatrixROmem{self.M==N,self.N==N} {     
        super(N, N); // DO NOT PUT ANYTHING BEFORE THIS LINE!!!
        Console.OUT.printf("GMatrixROmem.x10 initialization\n");
        this.bfs = bfs;
        this.mol = molecule;
        this.nOrbital = nOrbital;

        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        if (jd.roNK==-1) this.roNK=jd.roN;
        else this.roNK=jd.roNK;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
  
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N, N);
        kMatrix = new DenseMatrix(N, N);
        scratch = new DenseMatrix(N, N);

        val roLm = (roL+1)*(roL+1);
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roLm);
        aux = new Integral_Pack(roN,roL);
        dk = new Rail[Double](roLm); // eqn 15b in RO#7
        muk = new DenseMatrix(N,roLm); // Biggest RO array stored in Memory

        // Shell/Shellpair counting & allocation 
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }


        // Reconstruct Table IV in LHG2012 paper
        Console.OUT.printf("maxam=%d\n",maxam);
        val F1 = new Array[Int](0..(2*maxam));
        val F2 = new Array[Int](0..(2*maxam));
        val F3 = new Array[Int](0..maxam*0..maxam);
        val F4 = new Array[Int](0..maxam*0..maxam);
        val F2e = new Array[Int](0..(2*maxam));
        for (var a:Int=1; a<=2*maxam; a++) {
            F2e(a)=0;
            for (var i:Int=0; i<=a; i++) 
                for (var j:Int=0; j<=a-i; j++) {
                    val k=a-i-j;
                    if (k==1) F2e(a)+=3;
                    else if (j==1) F2e(a)+=4;
                    else if (i==1) F2e(a)+=4;
                    else if (k>=2) F2e(a)+=5;
                    else /* */ F2e(a)+=6;
                }
            Console.OUT.printf("F2e(%d)=%d\n",a,F2e(a));
        }        
        for (var a:Int=0; a<=maxam; a++) for (var b:Int=0; b<=a; b++) {
            F1(a+b)=3*(a+b+1);
            F2(a+b)=0;
            for (var f:Int=1; f<=a+b; f++)
                F2(a+b)+=(a+b-f+1)*F2e(f);
            F3(a,b)=nCr(a+b+3,3)-nCr(a+2,3);
            F4(a,b)=0;
            for (var f:Int=1; f<=b; f++) for (var e:Int=a; e<=a+b-f; e++)  
                F4(a,b)+=nCr(e+2,2)*nCr(f+2,2);
            F4(a,b)*=2;
            Console.OUT.printf("%2d %2d %5d %5d %5d %5d\n",a,b,F1(a+b),F2(a+b),F3(a,b),F4(a,b));          
        }

        // The cost factor here will be used later

        val threshold = 1.0e-8;
        val rawShellPairs = new Array[ShellPair](nShell*(nShell+1)/2); 
        val zeroPoint = Point3d(0.0, 0.0, 0.0);
        val emptyRailD = new Rail[Double](0),emptyRailI = new Rail[Int](0); 
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRailD, emptyRailD, emptyRailD, emptyRailD, 0, 0, 0, 0, emptyRailI, threshold);

        Console.OUT.printf("roN=%d roL=%d roNK=%d\n",roN,roL,roNK); 
        var mu:Int=0,nu:Int=0,ind:Int=0; 
        // centre a       
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
                val iaFunc = aFunc.get(i);
                // centre b
                for(var b:Int=0; b<noOfAtoms; b++) {
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = bFunc.size();
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);
                        //Console.OUT.printf("a=%d i=%d b=%d j=%d [naFunc=%d nbFunc=%d]\n", a,i,b,j,naFunc,nbFunc);
                        var aaFunc:ContractedGaussian=iaFunc,bbFunc:ContractedGaussian=jbFunc;
                        val aa=iaFunc.getTotalAngularMomentum();
                        val bb=jbFunc.getTotalAngularMomentum();
                        val maxbraa = (aa+1)*(aa+2)/2; 
                        val maxbrab = (bb+1)*(bb+2)/2;                  
                        if (aa>bb || (aa==bb && mu>=nu)) {                        
                            // extract info from basisfunctions
                            // Note that iaFunc and jbFunc are ContractedGaussians
                            // val may not work?  must specify the same type as in cpp code?
                            val aang = aaFunc.getTotalAngularMomentum();
                            val aPoint = aaFunc.origin;
                            val zetaA = aaFunc.exponents;
                            val conA = aaFunc.coefficients;
                            val dConA = conA.size;
                            val bang = bbFunc.getTotalAngularMomentum();
                            val bPoint = bbFunc.origin; 
                            val zetaB = bbFunc.exponents;                        
                            val conB = bbFunc.coefficients; 
                            val dConB = conB.size;
                            var contrib : Double = 0.; // S = conservative estimate
                            val R2 = Math.pow(aPoint.i-bPoint.i,2.)+Math.pow(aPoint.j-bPoint.j,2.)+Math.pow(aPoint.k-bPoint.k,2.);
                            for (var ii:Int=0; ii<dConA; ii++) for (var jj:Int=0; jj<dConB; jj++) 
                                contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2); // norm already included in con coef
                            val maxL = new Rail[Int](roN);
                            if (Math.abs(contrib)>threshold) 
                                rawShellPairs(ind++) = new ShellPair(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, mu, nu, maxL, Math.abs(contrib));                            
                        }
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab;
                        else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   
        numSigShellPairs = ind;
        shellPairs = new Array[ShellPair](numSigShellPairs); 
        for (i in 0..(numSigShellPairs-1))
            shellPairs(i)=rawShellPairs(i);

        val compareShellPairContribs = (y:ShellPair,x:ShellPair) => x.contrib.compareTo(y.contrib);
        ArrayUtils.sort[ShellPair](shellPairs, compareShellPairContribs);

        // numSigShellPairs = Math.abs(ArrayUtils.binarySearch[ShellPair](rawShellPairs, dummySignificantPair, compareShellPairContribs));
        Console.OUT.println("found " +/*+ numSigShellPairs + " significant from " + raw*/shellPairs.size /*+ " total"*/);

//      shellPairs = new Array[ShellPair](numSigShellPairs);
        var fhCost:Double=0.;
        var pAuxCount:Double=0.;
        var AuxCount:Double=0.;
        for (i in 0..(numSigShellPairs-1)) {

            val sh = /*raw*/shellPairs(i);
            val a=sh.aang; val b=sh.bang; val ka=sh.dconA; val kb=sh.dconB;

            // Find MaxL(n)
            var maxl:Int=0,maxn:Int=0,maxmaxl:Int=0; val THRESH=1.0e-8;                        
            var count1:Int=0; //sh.maxL = new Rail[Int](roN);
            for (var ron:Int=0; ron<=roN; ron++) {                           
                aux.genClass(sh.aang, sh.bang, sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.conA, sh.conB, sh.dconA, sh.dconB, temp, ron, roL);  

                var auxint:Double=0.,rol:Int=roL;                               
                for (; rol>=0 && auxint<THRESH; rol--) { 
                    for (var rom:Int=-rol; rom<=rol; rom++)  for (var tmu:Int=0; tmu<sh.maxbraa; tmu++) for (var tnu:Int=0; tnu<sh.maxbrab; tnu++) {
                        val ml=rol*(rol+1)+rom;
                        val mnu=tnu*(roL+1)*(roL+1)+ml;
                        val mmu=tmu*(roL+1)*(roL+1)*sh.maxbrab+mnu;
                        auxint = Math.max(Math.abs(temp(mmu)), auxint);
                    }
                }
                if (auxint<THRESH) maxl=-1; 
                else {
                    maxl=rol+1;
                    count1+=Math.pow(maxl+1,2);
                    maxn=ron;
                    if (maxl>maxmaxl) maxmaxl=maxl;
                }
                sh.maxL(ron)=maxl;
            }
            Console.OUT.printf("L(n=%d)=%d\n",ron,maxl);

            val c=(F1(a+b)+F2(a+b)+F3(a,b))*ka*kb+F4(a,b); 
            var K:Double=0;
            for (var ron:Int=0; ron<roN; ron++)
               K+=Math.pow(sh.maxL(ron)+1,2);
            fhCost+=c*K;
            val bracount=(a+1)*(a+2)/2*(b+1)*(b+2)/2;
            pAuxCount+=ka*kb*bracount*K;
            AuxCount+=bracount*K; 
        }
        Console.OUT.printf("nShell=%d nShellPairs=%d nSigShellPairs=%d fhCost=%e pAuxCount=%e AuxCount=%e\n",nShell,ind,numSigShellPairs,fhCost,pAuxCount,AuxCount);
    }

    private def nCr(a:Int,b:Int):Int {
        //return 1; 
        var re:Int=1;
        for (var temp:Int=b+1; temp<=a; temp++)
            re*=temp; // this can overflow easily - be careful
        for (var temp:Int=2; temp<=a-b; temp++)
            re/=temp; // this is always an integer - don't worry
        return re;
    }

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        Console.OUT.printf("GMatrixROmem.x10 compute\n");
        timer.start(0);
        jMatrix.reset();
        kMatrix.reset();

        // Form J matrix
        timer.start(1);
        var fac:Double; var ind:Int;
        for (var ron:Int=0; ron<=roN; ron++) {
            // Form dk vector
            dk.clear();
       
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (maxLron>=0) {
                    ind=0; fac=1.; if (sp.mu!=sp.nu) fac=2.0;
                    aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron);      
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++)  
                    for (var rol:Int=0; rol<=maxLron ; rol++) for (var rom:Int=-rol; rom<=rol; rom++)
                        dk(rol*(rol+1)+rom) += density(tmu,tnu)*fac*norm(tmu)*norm(tnu)*temp(ind++); // eqn 15b // skip some k's                 
                }
            }
             
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (sp.maxL(ron)>=0) { 
                    ind=0; 
                    aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron); 
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                        var jContrib:Double = 0.0;  
                        for (var rol:Int=0; rol<=maxLron ; rol++) for (var rom:Int=-rol; rom<=rol; rom++)
                            jContrib+=dk(rol*(rol+1)+rom)*norm(tmu)*norm(tnu)*temp(ind++);
                        jMatrix(tmu,tnu) += jContrib;   
                        if (sp.mu!=sp.nu) jMatrix(tnu,tmu) += jContrib;  
                    } 
                }
            }            
        }     

        timer.stop(1);
// vvvv For development purpose vvvvvv
        val eJ = scratch.mult(density, jMatrix).trace();
        Console.OUT.printf("  EJ = %.6f a.u.\n", eJ);
// ^^^^ It is not required for normal calculation ^^^^^
        Console.OUT.printf("    Time to construct JMatrix with RO: %.3g seconds\n", (timer.last(1) as Double) / 1e9);

        // Form K matrix
        timer.start(2);

        if (counter++!=0) // First cycle gives EK=0 
        for(aorb in 0..(nOrbital-1)) for (var ron:Int=0; ron<=roNK; ron++){ // To save mem
            muk.reset();
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (maxLron>=0) {
                    aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron); 
                    ind=0;                   
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                        val normMo = mos(aorb,tnu) * norm(tmu) * norm(tnu);
                        for (var rol:Int=0; rol<=maxLron; rol++) for (var rom:Int=-rol; rom<=rol; rom++)
                            muk(tmu,rol*(rol+1)+rom) += normMo * temp(ind++); // skip some k's    
                    }                                   
                    if (sp.mu!=sp.nu) {
                        ind=0;
                        for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                            val normMo = mos(aorb,tmu)*norm(tmu)*norm(tnu);
                            for (var rol:Int=0; rol<=maxLron; rol++) for (var rom:Int=-rol; rom<=rol; rom++) 
                                muk(tnu,rol*(rol+1)+rom) += normMo * temp(ind++); // skip some k's    
                        }
                    } 
                }              
            }

            kMatrix.multTrans(muk, muk, true);
        }   
        
        timer.stop(2);
// vvvv For development purpose vvvvvv
        val eK = scratch.mult(density, kMatrix).trace();
        Console.OUT.printf("  EK = %.6f a.u.\n", eK);
// ^^^^ It is not required for normal calculation ^^^^^
        Console.OUT.printf("    Time to construct KMatrix with RO: %.3g seconds\n", (timer.last(2) as Double) / 1e9);

        // Form G matrix
        jMatrix.d.map(this.d, kMatrix.d, (j:Double,k:Double)=>(2.0*j-k)); // eqn 14
        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(0) as Double) / 1e9);

    }
 
}


