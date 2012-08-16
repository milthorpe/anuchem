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
import x10.matrix.blas.DenseMatrixBLAS;

import x10x.vector.Vector;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.qm.ShellPair; 
import au.edu.anu.qm.Ylm; 
import au.edu.anu.qm.AuxInt; 
import au.edu.anu.qm.IntLm; 
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;
import edu.utk.cs.papi.PAPI;

public class GMatrixROmem2 extends DenseMatrix{self.M==self.N} {
    public val timer = new StatisticalTimer(5);
    static TIMER_TOTAL = 0;
    static TIMER_JMATRIX = 1;
    static TIMER_KMATRIX = 2;
    static TIMER_GENCLASS = 3;

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];
  
    val nOrbital:Int;
    val numSigShellPairs:Int;
    val omega:Double;    
    val norm:Rail[Double];
    val temp:Rail[Double];
    val dk:Rail[Double];
    val emptyRailD = new Rail[Double](0),emptyRailI = new Rail[Int](0); 

    val roN:Int; val roNK:Int;  val roL:Int; val roZ:Double;

    // Big arrays stored in Memory
    val shellPairs:Rail[ShellPair]; 
    val auxInts:Rail[AuxInt]; 
    val muk:DenseMatrix{self.M==this.N}; 
    val maxmaxl:Rail[Int];
    val jMatrix:DenseMatrix{self.M==self.N,self.N==this.N};
    val kMatrix:DenseMatrix{self.M==self.N,self.N==this.N};

    // used for calculating eJ, eK
    val scratch:DenseMatrix{self.M==self.N,self.N==this.N};

    // PAPI performance counters
    transient val papi:PAPI;

    transient val aux:Integral_Pack;
    var counter:Int=0;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int, omega:Double, roThresh:Double):GMatrixROmem2{self.M==N,self.N==N} {     
        super(N, N); 
        val result = Runtime.execForRead("date"); Console.OUT.printf("GMatrixROmem.x10 initialization starts at %s...\n",result.readLine()); 
        this.bfs = bfs;
        this.mol = molecule;
        this.nOrbital = nOrbital;
        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        if (jd.roNK==-1) this.roNK=jd.roN; else this.roNK=jd.roNK;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
        this.omega = jd.omega;
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N, N);
        kMatrix = new DenseMatrix(N, N);
        scratch = new DenseMatrix(N, N);
        val roLm = (roL+1)*(roL+1);
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roLm); // will be passed to C++ code
        aux = new Integral_Pack(roN,roL,omega);
        dk = new Rail[Double](roLm); // eqn 15b in RO#7
        muk = new DenseMatrix(N,roLm); 

        // Reconstruct Table IV in LHG2012 paper
        Console.OUT.printf("maxam=%d\n",maxam);
        val F1 = new Array[Int](0..(2*maxam)); val W1 = new Array[Int](0..(2*maxam));
        val F2 = new Array[Int](0..(2*maxam)); val W2 = new Array[Int](0..(2*maxam));
        val F3 = new Array[Int](0..maxam*0..maxam); val W3 = new Array[Int](0..maxam*0..maxam);
        val F4 = new Array[Int](0..maxam*0..maxam); val W4 = new Array[Int](0..maxam*0..maxam);
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
                    F2e(a)++; // REALITY FACTOR
                }
            // Console.OUT.printf("F2e(%d)=%d\n",a,F2e(a));
        }      
        Console.OUT.printf(" a  b  F123    F4  W123    W4\n-------------------------------\n");  
        for (var a:Int=0; a<=maxam; a++) for (var b:Int=0; b<=a; b++) {
            F1(a+b)=3*(a+b+1); W1(a+b)=a+b+1;
            F2(a+b)=0; W2(a+b)=nCr(a+b+4,4);
            for (var f:Int=1; f<=a+b; f++)
                F2(a+b)+=(a+b-f+1)*F2e(f);
            if (a>0) F3(a,b)=nCr(a+b+3,3)-nCr(a+2,3); else F3(0,0)=1;
            F3(a,b)*=2; // REALITY FACTOR
            W3(a,b)=F3(a,b);
            F4(a,b)=0;
            for (var f:Int=1; f<=b; f++) for (var e:Int=a; e<=a+b-f; e++)  
                F4(a,b)+=nCr(e+2,2)*nCr(f+2,2);
            W4(a,b)=F4(a,b); F4(a,b)*=2;
            Console.OUT.printf("%2d %2d %5d %5d %5d %5d\n",a,b,F1(a+b)+F2(a+b)+F3(a,b),F4(a,b),W1(a+b)+W2(a+b)+W3(a,b),W4(a,b));          
            //Console.OUT.printf("%2d %2d %5d %5d %5d %5d\n",a,b,F1(a+b),F2(a+b),F3(a,b),F4(a,b));  
        }

        // Shell/Shellpair business 
        Console.OUT.printf("roN=%d roL=%d roNK=%d\n",roN,roL,roNK); 
        val threshold = 1.0e-8; // ******************
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }
        val rawShellPairs = new Array[ShellPair](nShell*(nShell+1)/2); 
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
                                contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*R2); 
                            contrib=Math.abs(contrib);
                            if (contrib>=threshold) {
                                val maxL = new Rail[Int](roN);
                                rawShellPairs(ind++) = new 
                                    ShellPair(aang,bang,aPoint,bPoint,zetaA,zetaB,conA,conB,dConA,dConB,mu,nu,maxL,contrib);             
                            }               
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
            shellPairs(i)=rawShellPairs(i); // TODO delete rawShellPairs
        // val compareShellPairContribs = (y:ShellPair,x:ShellPair) => x.contrib.compareTo(y.contrib);
        // ArrayUtils.sort[ShellPair](shellPairs, compareShellPairContribs);
        // numSigShellPairs = Math.abs(ArrayUtils.binarySearch[ShellPair](rawShellPairs, dummySignificantPair, compareShellPairContribs));

        // Find MaxL(n) -- This screening is a gold standard but slow
                               
        var fCost:Double=0.;
        var wCost:Double=0.;
        var mCost:Double=0.;
        var pAuxCount:Double=0.;
        var AuxCount:Double=0.;
   
        auxInts = new Rail[AuxInt](numSigShellPairs);
        for (i in 0..(numSigShellPairs-1)) {
            val sh = shellPairs(i); val a=sh.aang; val b=sh.bang; val ka=sh.dconA; val kb=sh.dconB;            
            var maxl:Int=0,maxn:Int=0; 
            var count1:Int=0; 
            val intLms = new Rail[IntLm](roN); 
            val tempY = new Rail[Double](sh.dconA*sh.dconB*(roL+1)*(roL+1));
            aux.genClassY(sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.dconA, sh.dconB,  roL, tempY);

            for (var ron:Int=0; ron<=roN; ron++) {                           
                aux.genClass(sh.aang, sh.bang, sh.aPoint, sh.bPoint, sh.zetaA, sh.zetaB, sh.conA, sh.conB, sh.dconA, sh.dconB, temp, ron, roL, tempY, roL);  
                var auxint:Double=-1/*can't set to 0 as roThresh can be 0*/,rol:Int=roL;                               
                for (; rol>=0 && auxint<roThresh; rol--) { 
                    for (var rom:Int=-rol; rom<=rol; rom++)  for (var tmu:Int=0; tmu<sh.maxbraa; tmu++) for (var tnu:Int=0; tnu<sh.maxbrab; tnu++) {
                        val ml=rol*(rol+1)+rom;
                        val mnu=tnu*(roL+1)*(roL+1)+ml;
                        val mmu=tmu*(roL+1)*(roL+1)*sh.maxbrab+mnu;
                        auxint = Math.max(Math.abs(temp(mmu)), auxint);
                    }
                }
                if (auxint<roThresh) maxl=-1; 
                else {
                    maxl=rol+1;
                    count1+=Math.pow(maxl+1,2);
                    maxn=ron;                   
                }
                sh.maxL(ron)=maxl; 
                if (maxl>=0) {
                    // copy temp to arr
                    val arr = new Rail[Double]((maxl+1)*(maxl+1)*sh.maxbraa*sh.maxbrab); mCost+=(maxl+1)*(maxl+1)*sh.maxbraa*sh.maxbrab;
                    ind=0;
                    for (var tmu:Int=0; tmu<sh.maxbraa; tmu++) for (var tnu:Int=0; tnu<sh.maxbrab; tnu++) 
                        for (rol=0; rol<=maxl ; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                            val ml=rol*(rol+1)+rom;
                            val mnu=tnu*(roL+1)*(roL+1)+ml;
                            val mmu=tmu*(roL+1)*(roL+1)*sh.maxbrab+mnu;
                            arr(ind++)=temp(mmu);
                        }
                    intLms(ron) = new IntLm(arr,maxl); 
                }
                else intLms(ron) = new IntLm(emptyRailD,-1);
            }
            auxInts(i) = new AuxInt(intLms,maxn);
            val fc=(F1(a+b)+F2(a+b)+F3(a,b))*ka*kb+F4(a,b); 
            val wc=(F1(a+b)/3.+nCr(a+b+4,4)+F3(a,b))*ka*kb+F4(a,b)/2.; 
            var K:Double=0;
            for (var ron:Int=0; ron<roN; ron++)
               if (sh.maxL(ron)>=0) K+=Math.pow(sh.maxL(ron)+1,2);
            fCost+=fc*K;
            wCost+=wc*K;
            val bracount=(a+1)*(a+2)/2*(b+1)*(b+2)/2;
            pAuxCount+=ka*kb*bracount*K;
            AuxCount+=bracount*K; 
        }
        maxmaxl = new Rail[Int](roN);
        for (var ron:Int=0; ron<=roN; ron++) {
            maxmaxl(ron)=-1;
            for (i in 0..(numSigShellPairs-1))
                if (shellPairs(i).maxL(ron)>maxmaxl(ron)) maxmaxl(ron)=shellPairs(i).maxL(ron);
            Console.OUT.printf("maxmaxl(%d)=%d\n", ron, maxmaxl(ron));
        }

        Console.OUT.printf("nShell=%d numSigShellPairs=%d pAuxCount=%e (primitive) AuxCount=%e (contracted)\n",nShell,numSigShellPairs,pAuxCount,AuxCount);
        Console.OUT.printf("mCost=%e\n", mCost);
        Console.OUT.printf("fCost=%e\n", fCost);
        Console.OUT.printf("wCost=%e\n", wCost);
        papi = new PAPI();
        //papi.countFlops();
        papi.countMemoryOps();
    }

    private def nCr(a:Int,b:Int):Int {
        var re:Int=1;
        for (var t:Int=b+1; t<=a; t++)
            re*=t; // this can overflow easily - be careful
        for (var t:Int=2; t<=a-b; t++)
            re/=t; // this is always an integer - don't worry
        return re;
    }

    public def compute(density:Density{self.N==this.N}, mos:MolecularOrbitals{self.N==this.N}) {
        val result = Runtime.execForRead("date"); 
        Console.OUT.printf("GMatrixROmem.x10 compute starts at %s...\n",result.readLine()); 
        timer.start(TIMER_TOTAL);
        jMatrix.reset();
        kMatrix.reset();

        var prevLoadCount:Long = 0;

        // Form J matrix
        papi.reset();
        timer.start(TIMER_JMATRIX); var t:Double=0.;
        var fac:Double; var ind:Int; var jContrib:Double;

        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) fac=2.; else fac=1.;
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++)
                scratch(tmu,tnu)=density(tmu,tnu)*fac*norm(tmu)*norm(tnu);
        }
        if (counter==0)  Console.OUT.printf("rawtime\n");
        for (var ron:Int=0; ron<=roN; ron++) {
            // Form dk vector
            dk.clear();       
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);                
                if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1);
                    ind=0;  
                    timer.start(TIMER_GENCLASS);
                    val temp1=auxInts(spInd).intLms(ron).IntLm;
                    //papi.start();
                    //aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron,ylms(spInd).y, ylms(spInd).maxL);
                    //papi.stop();
                    timer.stop(TIMER_GENCLASS);
                    if (counter==0) {
                        //Console.OUT.printf("%d\t%d\t%d\t%d\t%d\t%d",sp.aang, sp.bang, sp.dconA, sp.dconB, ron, maxLron); // printf can accomodate upto 6 arguments?
                        //Console.OUT.printf("\t%20.15e\n",(timer.last(TIMER_GENCLASS) as Double)/1e9);
                        /*val loadCount = papi.getCounter(1);
                        Console.OUT.printf("\t%lld\n",loadCount-prevLoadCount);
                        prevLoadCount = loadCount;*/
                    }

                    t+=(timer.last(TIMER_GENCLASS) as Double)/1e9;
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                        val scdmn=scratch(tmu,tnu);
                        for (var rolm:Int=0; rolm<maxLm; rolm++)
                            dk(rolm) += scdmn*temp1(ind++); // eqn 15b 
                    }
                }
            }
             
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (sp.maxL(ron)>=0) { 
                    val maxLm=(maxLron+1)*(maxLron+1);
                    ind=0;  
                    //timer.start(TIMER_GENCLASS);
                    //papi.start();
                    val temp1=auxInts(spInd).intLms(ron).IntLm;
                    //aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron,ylms(spInd).y, ylms(spInd).maxL);  
                    //papi.stop();
                    //timer.stop(TIMER_GENCLASS); t+= (timer.last(TIMER_GENCLASS) as Double)/1e9 ;

                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                        val nrm=norm(tmu)*norm(tnu);
                        jContrib = 0.;
                        for (var rolm:Int=0; rolm<maxLm; rolm++)
                            jContrib+=dk(rolm)*nrm*temp1(ind++);
                        jMatrix(tmu,tnu) += jContrib;
                    } 

                }
            }            
        } 
        if (counter==0) Console.OUT.printf("rawtime\n");    
        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            if (sp.mu!=sp.nu) for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) 
                jMatrix(tnu,tmu) = jMatrix(tmu,tnu);
        }
        timer.stop(TIMER_JMATRIX);
        Console.OUT.print("J matrix ");
        //papi.printFlops();
        papi.printMemoryOps();


// vvvv For development purpose vvvvvv
        val eJ = scratch.mult(density, jMatrix).trace();
        Console.OUT.printf("  EJ = %.6f a.u.\n", eJ);
// ^^^^ It is not required for normal calculation ^^^^^

        Console.OUT.printf("    Time to construct JMatrix with RO: %.3g seconds (%.4g for ints)\n", (timer.last(TIMER_JMATRIX) as Double) / 1e9, t);

        // Form K matrix
        // if (counter++!=0) // if guess MOs is not available - i.e. SAD, first cycle gives EK=0 // This is partly fixed by using 'core' guess for MOs
        //papi.resetCount();
        Console.OUT.print("K matrix ");
        timer.start(TIMER_KMATRIX); t=0.;
        var t1:Double=0.;  var t2:Double=0.; var t3:Double=0.;
        // Version 3 - Unlimited memory - try to see what performance we can get
        val roK = (roL+1)*(roL+1);
        val auxIntMat = new DenseMatrix(N,N*roK);
        val halfAuxMat = new DenseMatrix(nOrbital,N*roK);

for (var ron:Int=0; ron<=roNK; ron++){

        // step 1 - could be moved to initialization // need more mem
        timer.start(TIMER_GENCLASS);
        
       
        for (spInd in 0..(numSigShellPairs-1)) {
            val sp=shellPairs(spInd);
            val maxLron=sp.maxL(ron);
            if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1); 
                    val temp1=auxInts(spInd).intLms(ron).IntLm; 
                    ind=0; 
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                        val normMuNu = norm(tmu)*norm(tnu);          
                        val baseMu = tmu*roK; val baseNu = tnu*roK;
                        for (var rolm:Int=0; rolm<maxLm; rolm++) {
                            auxIntMat(tmu,baseNu+rolm) = auxIntMat(tnu,baseMu+rolm) = normMuNu*temp1(ind);
                            ind++;
                        }
                    }
            }
        }
        timer.stop(TIMER_GENCLASS); t1+= (timer.last(TIMER_GENCLASS) as Double)/1e9 ;
       
        // step 2
        timer.start(TIMER_GENCLASS);
        
        // Matrix Matrix multiplication
        for (var aorb:Int=0;  aorb<nOrbital; aorb++) for (var muk:Int=0;  muk<N*roK; muk++) {
            halfAuxMat(aorb,muk)=0.;
            for (var lambda:Int=0;  lambda<N; lambda++)
                halfAuxMat(aorb,muk)+=mos(aorb,lambda)*auxIntMat(lambda,muk);
        }
        timer.stop(TIMER_GENCLASS); t2+= (timer.last(TIMER_GENCLASS) as Double)/1e9 ;

        // step 3
        timer.start(TIMER_GENCLASS);
        val halfAuxMat2 = new DenseMatrix(N,nOrbital*roK);
        for (var mu:Int=0;  mu<N; mu++)  for (var aorb:Int=0;  aorb<nOrbital; aorb++) for (var k:Int=0;  k<roK; k++) 
            halfAuxMat2(mu,aorb*roK+k)=halfAuxMat(aorb,mu*roK+k);
            
        // Matrix Matrix multiplication
        kMatrix.multTrans(halfAuxMat2, halfAuxMat2, true);
        timer.stop(TIMER_GENCLASS); t3+= (timer.last(TIMER_GENCLASS) as Double)/1e9 ;
        //for (var mu:Int=0;  mu<N; mu++) for (var nu:Int=0; nu<N; nu++) {
        //}
        
}
        Console.OUT.printf("version 3: time %lf %lf %lf \n",t1,t2,t3);
        // Version 2 - save memory and try to be smart by reducing numbers of loads
        /* 
        val giant3DMat = new DenseMatrix(N,(roL+1)*(roL+1)*nOrbital);
        //val t2DMat = new Array[Double](0..(nOrbital-1)*0..((roL+1)*(roL+1)-1));
        // prepare vector of MOs if LAPACK is to be used.
        for (var ron:Int=0; ron<=roNK; ron++){
            timer.start(TIMER_GENCLASS);
            giant3DMat.reset();
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1); 
                    val temp1=auxInts(spInd).intLms(ron).IntLm; 
                    ind=0; var ind2:Int=0;
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                        val normMo = norm(tmu)*norm(tnu);                        
                        // Vector Vector multiplication
                        ind2=0;
                        for (var rolm:Int=0; rolm<maxLm; rolm++) {
                            val normInt = normMo*temp1(ind++);
                            for (var aorb:Int=0;  aorb<nOrbital; aorb++) {
                                //t2DMat(aorb,rolm)=normInt*mos(aorb,tmu);
                                if (sp.mu!=sp.nu) {
                                    giant3DMat(tnu,ind2)+=normInt*mos(aorb,tmu);  // DAXPY incY = mos.N
                                }
                                giant3DMat(tmu,ind2++)+=normInt*mos(aorb,tnu);  // DAXPY incY = mos.N
                            }
                        }
                    }
                }
            }
            timer.stop(TIMER_GENCLASS); t+= (timer.last(TIMER_GENCLASS) as Double)/1e9;
            kMatrix.multTrans(giant3DMat, giant3DMat, true);
        } */
        // Version 1 - save as much as memory
        /* for (aorb in 0..(nOrbital-1)) for (var ron:Int=0; ron<=roNK; ron++){ // To save mem
            muk.reset();
            for (spInd in 0..(numSigShellPairs-1)) {
                val sp=shellPairs(spInd);
                val maxLron=sp.maxL(ron);
                if (maxLron>=0) {
                    val maxLm=(maxLron+1)*(maxLron+1);  
                    timer.start(TIMER_GENCLASS);
                    //papi.start();
                    // aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, ron, maxLron,ylms(spInd).y, ylms(spInd).maxL);
                    val temp1=auxInts(spInd).intLms(ron).IntLm;
                    //papi.stop();
                    
                    ind=0;                   
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                        val normMo = mos(aorb,tnu)*norm(tmu)*norm(tnu); 
                        for (var rolm:Int=0; rolm<maxLm; rolm++) 
                            muk(tmu,rolm)+= normMo*temp1(ind++);     
                    }                                   
                    if (sp.mu!=sp.nu) {
                        ind=0;
                        for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { 
                            val normMo = mos(aorb,tmu)*norm(tmu)*norm(tnu);
                             for (var rolm:Int=0; rolm<maxLm; rolm++) 
                                muk(tnu,rolm)+= normMo*temp1(ind++);    
                        }
                    } 
                    timer.stop(TIMER_GENCLASS); t+= (timer.last(TIMER_GENCLASS) as Double)/1e9;

                }              
            }        
            // kMatrix.multTrans(muk, muk, true); // most expensive -- can be reduced by using maxmaxl
            val m=(maxmaxl(roN)+1)*(maxmaxl(roN)+1);
            // val m = (roL+1)*(roL+1);
            val dim = [N, N, m];
            DenseMatrixBLAS.compMultTrans(muk, muk, kMatrix, dim, true);

        }   */
        timer.stop(TIMER_KMATRIX);
        
        //papi.printFlops();
        //papi.printMemoryOps();

// vvvv For development purpose vvvvvv
        val eK = scratch.mult(density, kMatrix).trace();
        Console.OUT.printf("  EK = %.6f a.u.\n", eK);
// ^^^^ It is not required for normal calculation ^^^^^

        Console.OUT.printf("    Time to construct KMatrix with RO: %.3g seconds (%.4g for ints)\n", (timer.last(TIMER_KMATRIX) as Double) / 1e9, t);

        // Form G matrix
        jMatrix.d.map(this.d, kMatrix.d, (j:Double,k:Double)=>(2.0*j-k)); // eqn 14
        timer.stop(TIMER_TOTAL);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(TIMER_TOTAL) as Double) / 1e9);

    }

}


