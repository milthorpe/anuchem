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

import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;

/**
 * Implementation of Hartree-Fock SCF method
 * @author: V. Ganesh, J. Milthorpe, T. Limpanuparb 
 */
public class HartreeFockSCFMethod extends SCFMethod { 

    val roZ:Double;

    public def this(mol:Molecule[QMAtom],  
                    oneE:OneElectronIntegrals, 
                    bfs:BasisFunctions) {
        super(mol, oneE, bfs);
        val jd = JobDefaults.getInstance();
        this.roZ=jd.roZ;
    }

    public def scf() : void {
        val timer = new StatisticalTimer(5);
        val TIMER_TOTAL:Int = 0; val TIMER_H:Int = 1; val TIMER_S:Int = 2; val TIMER_F:Int = 3; val TIMER_M:Int = 4; //Can't use static

        val noOfElectrons = molecule.getNumberOfElectrons();        
        if (noOfElectrons%2 != 0 || molecule.getMultiplicity()!=1) {
           Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        }
        val noOfOccupancies = noOfElectrons / 2;

        val nuclearEnergy = getNuclearEnergy();
        Console.OUT.printf("Nuclear repulsion energy = %.6f a.u.\n", nuclearEnergy);

        val noOfBasisFunctions:Long = bfs.getBasisFunctions().size();
        val noOfIntegrals:Long = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;
        Console.OUT.println("\nNumber of 2E integrals: " + noOfIntegrals);

        energy = 0.0;

        var converged:Boolean = false;
        var oldEnergy:Double = 0.0;
        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();
        @Ifdef("__DEBUG__") {
            val diag = new GMLDiagonalizer();
            diag.diagonalize(overlap);
            val eigval =  diag.getEigenValues().d;;
            Console.OUT.printf("Three smallest eigen values %e %e %e\n",eigval(0),eigval(1),eigval(2));        
            if (eigval(0)<1e-6) Console.OUT.printf("Linear dependence detected!!!\n");
        }
        // init memory for the matrices
        val N = hCore.N;
        val jd = JobDefaults.getInstance();

        val gMatrix:GMatrix{self.N==N};
        val gMatrixRo:GMatrixROmem2{self.N==N};
        val roThresh=jd.roThresh;
        val thresh=jd.thresh;
        if (jd.roOn>0 && maxIteration>0) {
            gMatrixRo = new GMatrixROmem2(N, bfs, molecule, noOfOccupancies,0.,roZ*roThresh);
        } else {
            gMatrixRo = null;
        }
        if ((jd.roOn==0 || jd.compareRo==true) && maxIteration>0) {
            gMatrix = new GMatrix(N, bfs, molecule,0.,roZ*thresh);
        } else {
            gMatrix = null;
        }

        val mos = new MolecularOrbitals(N);
        val density = new Density(N, noOfOccupancies); // density.make();

        var fock:Fock{self.M==N,self.N==N} = new Fock(N);

        if (jd.guess.equals(JobDefaults.GUESS_SAD)) {
            Console.OUT.printf("guess = SAD... (core for MOs)\n");
            density.applyGuess(bfs.getSAD());
            mos.compute(hCore, overlap); // Cheat
        } else if (jd.guess.equals(JobDefaults.GUESS_CORE)) {
            Console.OUT.printf("guess = core...\n");
            mos.compute(hCore, overlap);
            density.compute(mos);
        } else {
            Console.OUT.printf("guess = ???...\n");
        }
        Console.OUT.printf("Starting SCF procedure...\n");      

        val diis = new DIISFockExtrapolator(N);

        // start the SCF cycle
        for(var scfIteration:Int=0; scfIteration<maxIteration; scfIteration++) {
            // make or guess density
	    if (scfIteration>0) {
                density.compute(mos);
            }
            
            // make the G matrix
            if (jd.roOn>0) {
                gMatrixRo.compute(density, mos);           
                if (jd.compareRo==true) {
                    gMatrix.compute(density);
                    Console.OUT.println("G Mat");
                    Console.OUT.println(gMatrix);

                    Console.OUT.println("G Mat RO");
                    Console.OUT.println(gMatrixRo);
                }
            } else {
                gMatrix.compute(density);
            }            

            timer.start(TIMER_F);
            // make fock matrix
            if (jd.roOn>0) fock.compute(hCore, gMatrixRo);
            else fock.compute(hCore, gMatrix);
            // SCF_ALGORITHM = DIIS  
            fock = diis.next(fock, overlap, density);
            timer.stop(TIMER_F);
            Console.OUT.println ("    Time to construct Fock: " + (timer.total(TIMER_F) as Double) / 1e9 + " seconds");

            timer.start(TIMER_M);
            // compute the new MOs
            mos.compute(fock, overlap);
            timer.stop(TIMER_F);
            Console.OUT.println ("    Time to form MOS: " + (timer.total(TIMER_M) as Double) / 1e9 + " seconds");
         
            // compute the total energy 
            val eOne = density.clone().mult(density, hCore).trace();
            val eTwo = density.clone().mult(density, fock).trace();           
            energy = eOne + eTwo + nuclearEnergy;

            Console.OUT.printf("Cycle #%i Total energy = %.6f a.u. (scale factor = %.10f)", scfIteration, energy/roZ,roZ);
            if (scfIteration>0) {
                Console.OUT.printf(" (%.6f)", (energy-oldEnergy)/roZ);
            } else {
                // ignore the first cycle's timings 
                // as far fewer integrals are calculated
                if (jd.roOn==0 || jd.compareRo==true) {
                    gMatrix.timer.clear(0);
                }
                if (jd.roOn>0) {
                    gMatrixRo.timer.clear(0);
                }
            }
            Console.OUT.printf("\n----------------------------------------------------------\n");
            // check for convergence
            if (scfIteration > 0 && diis.isConverged()) {
                converged = true;
                break;
            } // end if
            
            oldEnergy = energy;
        } // end of SCF iteration

        if (!converged) 
           Console.OUT.println("SCF did not converge in " + maxIteration + " cycles!");
        else 
           Console.OUT.printf("SCF converged. Final SCF energy = %.6f a.u.\n", energy/roZ,roZ);

        Console.OUT.printf("==========================================================\n");

        var ecount:Double=0.; var maxDen:Double=0.; var minDen:Double=1e100;

        /* val pdmat = new DenseMatrix(N, N);
        pdmat.multTrans(density, overlap, true);
        for (var i:Int=0; i<N; i++) {
            ecount+=pdmat(i);
            if (pdmat(i)<0.) Console.OUT.printf("padmat(%d)=%e\n",i,pdmat(i));
        }      
        // Trace of PS = # electrons
        ecount=0.;*/
        for (var i:Int=0; i<N; i++) for (var j:Int=i; j<N; j++) {
             val contrib=density(i,j)*overlap(i,j);
             ecount+=contrib;
             if (i!=j) ecount+=contrib;
             if (density(i,j)>100.) Console.OUT.printf("density(%d,%d)=%e, s=%e\n",i,j,density(i,j),overlap(i,j));
             if (density(i,j)>maxDen) maxDen=density(i,j);
             else if (density(i,j)<minDen) minDen=density(i,j);
        }
        Console.OUT.printf("ecount=%e maxDen=%e minDen=%e\n",ecount,maxDen,minDen);

        if ((jd.roOn == 0 || jd.compareRo) && maxIteration>0) {
            Console.OUT.println("GMatrix construction timings:");
            gMatrix.timer.printSeconds();
        }

        if (jd.roOn>0 && maxIteration>0) {
            Console.OUT.println("GMatrixROmem construction timings:");
            gMatrixRo.timer.printSeconds();
        }
    
        // long range energy
        //Console.OUT.println("before RO heapSize = " + System.heapSize());
        while (jd.roOn>0 && jd.roN>0) {
            computeLongRangeRO(N, mos, noOfOccupancies, density, jd, bfs);
            System.gc();
            Console.OUT.print("Input new roN Omega roThresh (nnn ooo t) or 000 000 0 to exit:");
            val rbuf = new Rail[Byte](20);
            Console.IN.read(rbuf,0,10);
            jd.roN = ((rbuf(0)-48) as Int)*100+(rbuf(1)-48)*10+(rbuf(2)-48)*1;            
            jd.omega = (rbuf(4)-48)*1.0+(rbuf(5)-48)*0.1+(rbuf(6)-48)*0.01;
            jd.roThresh = Math.pow(10,-(rbuf(8)-48));
            Console.OUT.println("rbuf = "+rbuf);
            Console.OUT.println("new roN = "+jd.roN +" new omega ="+jd.omega+" new roThresh ="+jd.roThresh);
            
            //Console.OUT.println("after GC heapSize = " + System.heapSize());
        }

        if (jd.compareRo) {
            Console.OUT.println("Long-range - Conventional");
            val gMatrixL = new GMatrix(N, bfs, molecule,jd.roZ*jd.omega,roZ*jd.thresh); // RO Thesis Eq (2.22)
            gMatrixL.compute(density);   
        }
        Console.OUT.println("after conventional = " + System.heapSize());
            //fock.compute(hCore, gMatrixRo);
            //val eOne = density.clone().mult(density, hCore).trace();
            //val eTwo = density.clone().mult(density, fock).trace();
            //energy = eOne + eTwo + nuclearEnergy;
            //Console.OUT.printf("Cycle ** Total energy = %.6f a.u. (scale factor = %.6f)",  energy/roZ,roZ);
    }

    private def computeLongRangeRO(N:Int, mos:MolecularOrbitals{self.N==N}, noOfOccupancies:Int, density:Density{self.M==N,self.N==N}, jd:JobDefaults, bfs:BasisFunctions) {
         Console.OUT.println("Long-range - RO");
            density.compute(mos);
            val gMatrixRoL = new GMatrixROmem2(N, bfs, molecule, noOfOccupancies,jd.roZ*jd.omega,roZ*jd.roThresh); // RO Thesis Eq (2.22)
            gMatrixRoL.compute(density, mos);
/* 
            Console.OUT.println("after RO heapSize = " + System.heapSize() + " size of gMatrixRoL.muk = " + gMatrixRoL.muk.d.size
                 + " size of gMatrixRoL.auxIntMat = " + gMatrixRoL.auxIntMat.d.size
                 + " size of gMatrixRoL.jMatrix = " + gMatrixRoL.jMatrix.d.size);
*/
    }

}

