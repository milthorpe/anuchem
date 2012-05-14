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

import x10x.matrix.Matrix;
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

public class GMatrixROmem extends Matrix {

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    val roN:Int;
    val roL:Int;
    val roZ:Double;
    val nBasis:Int;
    val nOrbital:Int;
    
    val norm:Rail[Double];
    val temp:Rail[Double];
    val dk:Rail[Double];
    val muk:Array[Double](2){rect,zeroBased};

    val shellPairs:Array[ShellPair];

    transient val aux:Integral_Pack;
    val jMatrix:Matrix;
    val kMatrix:Matrix;

    var nSigShellPairs:Int;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], nOrbital:Int) {
        super(N);
        this.bfs = bfs;
        this.mol = molecule;
        this.nBasis=N;
        this.nOrbital = nOrbital;

        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
  
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new Matrix(N);
        kMatrix = new Matrix(N);

        val roK = (roN+1)*(roL+1)*(roL+1);
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roK);
        aux = new Integral_Pack(roN,roL);
        dk = new Rail[Double](roK); // eqn 15b in RO#7
        muk = new Array[Double](0..(nBasis-1)*0..(roK-1)); // Biggest RO array stored in Memory

        // Gen Shellpair
        var nShell:Int=0;
        val noOfAtoms = mol.getNumberOfAtoms();

        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            nShell+=aFunc.size();
        }

        shellPairs = new Array[ShellPair](nShell*(nShell+1)/2); 

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
                            val R = Math.sqrt(Math.pow(aPoint.i-bPoint.i,2.)+Math.pow(aPoint.j-bPoint.j,2.)+Math.pow(aPoint.k-bPoint.k,2.));
                            for (var ii:Int=0; ii<dConA; ii++) for (var jj:Int=0; jj<dConB; jj++) 
                                contrib+=conA(ii)*conB(jj)*Math.exp(-zetaA(ii)*zetaB(jj)/(zetaA(ii)+zetaB(jj))*Math.pow(R,2.)); // norm already included in con coef
                            // Console.OUT.printf("mu=%4d nu=%4d contrib=%17.10f\n",mu,nu,contrib);  

                            // Call genclass to find N and L appropriate to THRESH
                            aux.genClass(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, temp);
                            // Find Max N & Max L
                            var maxl:Int=0,maxn:Int=0; val THRESH=1.0e-7;var ron:Int; var rol:Int;
                           
                            var auxint:Double=0.;
                            for (var tmu:Int=0; tmu<maxbraa; tmu++) for (var tnu:Int=0; tnu<maxbrab; tnu++) {
                                for (ron=roN; ron>=0 && auxint<THRESH; ron--)                                    
                                    for (rol=0; rol<=roL; rol++) for (var rom:Int=-rol; rom<=rol; rom++) {
                                        val ml=rol*(rol+1)+rom;
                                        val mn=ron*(roL+1)*(roL+1)+ml;
                                        val mnu=tnu*(roN+1)*(roL+1)*(roL+1)+mn;
                                        val mmu=tmu*(roN+1)*(roL+1)*(roL+1)*maxbrab+mnu;
                                        if (Math.abs(temp(mmu))>auxint) auxint=Math.abs(temp(mmu));
                                    }
                                if (ron+1>maxn) maxn=ron+1; if (maxn==roN+1) maxn--;
                            }

                            auxint=0.;
                            for (var tmu:Int=0; tmu<maxbraa; tmu++) for (var tnu:Int=0; tnu<maxbrab; tnu++) {
                                for (rol=roL; rol>=0 && auxint<THRESH; rol--)                                    
                                    for (ron=0; ron<=roN; ron++) for (var rom:Int=-rol; rom<=rol; rom++) {
                                        val ml=rol*(rol+1)+rom;
                                        val mn=ron*(roL+1)*(roL+1)+ml;
                                        val mnu=tnu*(roN+1)*(roL+1)*(roL+1)+mn;
                                        val mmu=tmu*(roN+1)*(roL+1)*(roL+1)*maxbrab+mnu;
                                        if (Math.abs(temp(mmu))>auxint) auxint=Math.abs(temp(mmu));
                                    }
                                if (rol+1>maxl) maxl=rol+1;  if (maxl==roL+1) maxl--;
                            }
                            Console.OUT.printf("mu=%4d nu=%4d maxn=%4d maxl=%4d thresh=%e\n",mu,nu,maxn,maxl,THRESH);
                            shellPairs(ind++) = new ShellPair(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, mu, nu, roN, roL,Math.abs(contrib));
                        }
                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab;
                        else {mu+=maxbraa; nu=0;}
                    }    
                }
            }   
        }   

        new ArrayUtils[ShellPair]().sort(shellPairs, (y:ShellPair,x:ShellPair) => x.contrib.compareTo(y.contrib));
        for (nSigShellPairs=0; nSigShellPairs<=ind && shellPairs(nSigShellPairs).contrib > 1.0e-5; ) {
             val sh = shellPairs(nSigShellPairs);
             Console.OUT.printf("mu=%4d nu=%4d contrib=%17.10f\n",sh.mu,sh.nu,sh.contrib);
             nSigShellPairs++;
        }
        Console.OUT.printf("nShell=%d nShellPairs=%d nSigShellPairs=%d\n",nShell,ind,nSigShellPairs); 
    }


    public def compute(density:Density, mos:MolecularOrbitals) {
        timer.start(0);

        val N = getRowCount();
        val mosMat = mos.getMatrix();
        val denMat = density.getMatrix();
        val roK = (roN+1)*(roL+1)*(roL+1);

        // Form dk vector
        dk.clear();
        for (var spInd:Int=0; spInd<nSigShellPairs; spInd++) {
            val sp=shellPairs(spInd);
            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp);      
            var ind:Int=0; var fac:Double=1.;
            if (sp.mu!=sp.nu) fac=2.0;
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) for (var k:Int=0; k<roK; k++) 
               dk(k) += denMat(tmu,tnu)*fac*norm(tmu)*norm(tnu)*temp(ind++); // eqn 15b                   
        }     
    
        // Form J matrix
        val jMat = jMatrix.getMatrix();
        for (var spInd:Int=0; spInd<nSigShellPairs; spInd++) {
            val sp=shellPairs(spInd);
            aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp); 
    
            var ind:Int=0; 
            for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) {
                var jContrib:Double = 0.0;  
                for (var k:Int=0; k<roK; k++) 
                    jContrib+=dk(k)* norm(tmu)*norm(tnu)*temp(ind++);
                jMat(tmu,tnu) = jMat(tnu,tmu) = jContrib;   
            }            
        }     

        // Form K matrix
        val kMat = kMatrix.getMatrix();
        kMat.clear();

        for(aorb in 0..(nOrbital-1)) { // To save mem
            muk.clear();
            for (var spInd:Int=0; spInd<nSigShellPairs; spInd++) {
                val sp=shellPairs(spInd);
                aux.genClass(sp.aang, sp.bang, sp.aPoint, sp.bPoint, sp.zetaA, sp.zetaB, sp.conA, sp.conB, sp.dconA, sp.dconB, temp); 
                var ind:Int=0;          
                
                for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) for (var k:Int=0; k<roK; k++) 
                    muk(tmu,k) += mosMat(aorb,tnu) *norm(tmu)*norm(tnu)*temp(ind++);                                           
                if (sp.mu!=sp.nu) {
                ind=0;
                for (var tmu:Int=sp.mu; tmu<sp.mu+sp.maxbraa; tmu++) for (var tnu:Int=sp.nu; tnu<sp.nu+sp.maxbrab; tnu++) for (var k:Int=0; k<roK; k++) 
                    muk(tnu,k) += mosMat(aorb,tmu) *norm(tmu)*norm(tnu)*temp(ind++);
                }                
            }
            for (tmu in 0..(nBasis-1)) for (tnu in 0..(nBasis-1)) for (var k:Int=0; k<roK; k++) {
                //var kContrib:Double = 0.0;
                kMat(tmu,tnu) += muk(tmu,k)*muk(tnu,k); // check this statement
            }
        }   
        
        // Form G matrix
        val gMat = getMatrix();
        jMat.map(gMat, kMat, (j:Double,k:Double)=>(2.0*j-k)); // eqn 14

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }
 
}

