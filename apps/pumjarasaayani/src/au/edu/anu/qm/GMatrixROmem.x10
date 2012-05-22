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
import au.edu.anu.qm.ShellPair; // New
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;

/**
 * G matrix in HF calculation -- RO 
 */

public class GMatrixROmem extends DenseMatrix {

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    val roN:Int;
    val roL:Int;
    val roZ:Double;
    val nOrbital:Int;
    
    val norm:Rail[Double];
    val temp:Rail[Double];
    val dk:Rail[Double];
    val muk:DenseMatrix;

    val shellPairs:Rail[ShellPair];

    transient val aux:Integral_Pack;
    val jMatrix:DenseMatrix;
    val kMatrix:DenseMatrix;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int) {
        super(N, N);
        this.bfs = bfs;
        this.mol = molecule;
        this.nOrbital = nOrbital;

        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
  
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N, N);
        kMatrix = new DenseMatrix(N, N);

        val roK = (roN+1)*(roL+1)*(roL+1);
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roK);
        aux = new Integral_Pack(roN,roL);
        dk = new Rail[Double](roK); // eqn 15b in RO#7
        muk = new DenseMatrix(N, roK); // Biggest RO array stored in Memory

        // Gen Shellpair
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();

        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }

        val rawShellPairs = new Array[ShellPair](nShell*(nShell+1)/2); 

        var mu:Int = 0; 
        var nu:Int = 0; 
        // centre a
        var ind:Int=0;
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
                            // Console.OUT.printf("mu=%4d nu=%4d contrib=%17.10f\n",mu,nu,contrib);  

                            // Call genclass to find N and L appropriate to THRESH
                            aux.genClass(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, temp, roN, roL);
                            // Find Max N & Max L
                            var maxl:Int=0,maxn:Int=0; val THRESH=1.0e-10;var ron:Int; var rol:Int;
                           
                            var auxint:Double=0.; var count1:Int=0;
                            for (ron=roN; ron>=0 && auxint<THRESH; ron--)   {
                                for (var tmu:Int=0; tmu<maxbraa; tmu++) for (var tnu:Int=0; tnu<maxbrab; tnu++)                                  
                                    for (rol=0; rol<=roL; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                                        val ml=rol*(rol+1)+rom;
                                        val mn=ron*(roL+1)*(roL+1)+ml;
                                        val mnu=tnu*(roN+1)*(roL+1)*(roL+1)+mn;
                                        val mmu=tmu*(roN+1)*(roL+1)*(roL+1)*maxbrab+mnu;
                                        if (Math.abs(temp(mmu))>auxint) auxint=Math.abs(temp(mmu));
                                    }
                            }
                            maxn=ron+1; 
                            if (maxn==roN+1) { maxn--;  }

                            auxint=0.;
                            for (rol=roL; rol>=0 && auxint<THRESH; rol--)  {
                                for (var tmu:Int=0; tmu<maxbraa; tmu++) for (var tnu:Int=0; tnu<maxbrab; tnu++)                                  
                                    for (ron=0; ron<=roN; ron++) for (var rom:Int=-rol; rom<=rol; rom++) {
                                        val ml=rol*(rol+1)+rom;
                                        val mn=ron*(roL+1)*(roL+1)+ml;
                                        val mnu=tnu*(roN+1)*(roL+1)*(roL+1)+mn;
                                        val mmu=tmu*(roN+1)*(roL+1)*(roL+1)*maxbrab+mnu;
                                        if (Math.abs(temp(mmu))>auxint) auxint=Math.abs(temp(mmu)); 
                                    }
                            }
                            maxl=rol+1;  
                            if (maxl==roL+1) { maxl--;  }

                            Console.OUT.printf("mu=%4d nu=%4d THRESH=%e maxn=%d maxl=%d\n",mu,nu,THRESH,maxn,maxl);
                            rawShellPairs(ind++) = new ShellPair(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, mu, nu, maxn, maxl ,Math.abs(contrib));
                           
                            for (rol=0; rol<=roL; rol++) {
                                auxint=0.;
                                for (ron=roN; ron>=0 && auxint<THRESH; ron--)  {
                                    for (var rom:Int=-rol; rom<=rol; rom++) for (var tmu:Int=0; tmu<maxbraa; tmu++) for (var tnu:Int=0; tnu<maxbrab; tnu++) {
                                        val ml=rol*(rol+1)+rom;
                                        val mn=ron*(roL+1)*(roL+1)+ml;
                                        val mnu=tnu*(roN+1)*(roL+1)*(roL+1)+mn;
                                        val mmu=tmu*(roN+1)*(roL+1)*(roL+1)*maxbrab+mnu;
                                        if (Math.abs(temp(mmu))>auxint) auxint=Math.abs(temp(mmu));
                                    }
                                }
                                if (ron!=roN) maxn=ron+1; else maxn=roN; 
                                Console.OUT.printf("L=%d N=%d\n",rol,maxn);
                                count1+=(maxn+1)*Math.pow(rol+1,2);
                            }
                            Console.OUT.printf("K=%d\n",count1);


                            
                        }
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab;
                        else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   

        val compareShellPairContribs = (y:ShellPair,x:ShellPair) => x.contrib.compareTo(y.contrib);
        ArrayUtils.sort[ShellPair](rawShellPairs, compareShellPairContribs);
        val threshold = 1.0e-5;

        val zeroPoint = Point3d(0.0, 0.0, 0.0);
        val emptyRail = new Rail[Double](0);
        val dummySignificantPair = new ShellPair(0, 0, zeroPoint, zeroPoint, emptyRail, emptyRail, emptyRail, emptyRail, 0, 0, 0, 0, 0, 0, threshold);

        val numSigShellPairs = Math.abs(ArrayUtils.binarySearch[ShellPair](rawShellPairs, dummySignificantPair, compareShellPairContribs));
        Console.OUT.println("found " + numSigShellPairs + " significant from " + rawShellPairs.size + " total");

        shellPairs = new Array[ShellPair](numSigShellPairs);
        for (i in 0..(numSigShellPairs-1)) {
            val sh = rawShellPairs(i);
            Console.OUT.printf("mu=%4d nu=%4d contrib=%17.10f\n",sh.mu,sh.nu,sh.contrib);
            shellPairs(i) = sh;
        }
        Console.OUT.printf("nShell=%d nShellPairs=%d nSigShellPairs=%d\n",nShell,ind,numSigShellPairs);
    }


    public def compute(density:Density, mos:MolecularOrbitals) {
        timer.start(0);

        val roK = (roN+1)*(roL+1)*(roL+1);

        // Form dk vector
        dk.clear();
        for (spInd in 0..(shellPairs.size-1)) {
            val sp=shellPairs(spInd);
            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, sp.N, sp.L);      
            var ind:Int=0; var fac:Double=1.; var ron:Int; var rol:Int;
            if (sp.mu!=sp.nu) fac=2.0; 
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) //for (var k:Int=0; k<roK; k++) 
                for (ron=0; ron<=sp.N; ron++) for (rol=0; rol<=sp.L ; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                    val ml=rol*(rol+1)+rom;
                    val mn=ron*(roL+1)*(roL+1)+ml;
                    dk(mn) += density(tmu,tnu)*fac*norm(tmu)*norm(tnu)*temp(ind++); // eqn 15b // skip some k's                 
                }
        }     
    
        // Form J matrix
        for (spInd in 0..(shellPairs.size-1)) {
            val sp=shellPairs(spInd);
            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, sp.N, sp.L); 
    
            var ind:Int=0; var ron:Int; var rol:Int;
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                var jContrib:Double = 0.0;  
                // for (var k:Int=0; k<roK; k++) 
                //     jContrib+=dk(k)* norm(tmu)*norm(tnu)*temp(ind++); // skip some k's
                for (ron=0; ron<=sp.N; ron++) for (rol=0; rol<=sp.L ; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                    val ml=rol*(rol+1)+rom;
                    val mn=ron*(roL+1)*(roL+1)+ml;
                    jContrib+=dk(mn)* norm(tmu)*norm(tnu)*temp(ind++);
                }
                jMatrix(tmu,tnu) = jMatrix(tnu,tmu) = jContrib;   
            }            
        }     

        // Form K matrix
        kMatrix.reset();

        for(aorb in 0..(nOrbital-1)) { // To save mem
            muk.reset();
            for (spInd in 0..(shellPairs.size-1)) {
                val sp=shellPairs(spInd);
                aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp, sp.N, sp.L); 
                var ind:Int=0; var ron:Int; var rol:Int;         
                
                for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { //for (var k:Int=0; k<roK; k++)
                    val normMo = mos(aorb,tnu) * norm(tmu) * norm(tnu);
                    for (ron=0; ron<=sp.N; ron++) for (rol=0; rol<=sp.L ; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                        val ml=rol*(rol+1)+rom;
                        val mn=ron*(roL+1)*(roL+1)+ml;
                        muk(tmu,mn) += normMo * temp(ind++); // skip some k's    
                    }
                }                                   
                if (sp.mu!=sp.nu) {
                    ind=0;
                    for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) { //for (var k:Int=0; k<roK; k++) 
                        val normMo = mos(aorb,tmu) * norm(tmu) * norm(tnu);
                        for (ron=0; ron<=sp.N; ron++) for (rol=0; rol<=sp.L ; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                            val ml=rol*(rol+1)+rom;
                            val mn=ron*(roL+1)*(roL+1)+ml;
                            muk(tnu,mn) += normMo * temp(ind++); // skip some k's    
                        }                              
                       // muk(tnu,k) += mos(aorb,tmu) *norm(tmu)*norm(tnu)*temp(ind++); // skip some k's
                    }
                }                
            }

            kMatrix.multTrans(muk, muk, true);
        }   
        
        // Form G matrix
        jMatrix.d.map(this.d, kMatrix.d, (j:Double,k:Double)=>(2.0*j-k)); // eqn 14

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }
 
}

