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

import x10.io.File;
import x10.xrx.Runtime;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;

import au.edu.anu.qm.mta.Fragment; 
import au.edu.anu.qm.mta.Fragmentor; 

/**
 * PumjaRasaayani .. first, fully working, Quantum Chemistry code written in X10 ;-)
 *
 * Pumja : Sanskrit for Quantum
 * Rasaayn : Chemical (Rasayan Shastra : Chemistry / Chemical Science)
 *
 * @author: V.Ganesh
 */
public class PumjaRasaayani {
    var mol:Molecule[QMAtom];
    val inputFileName:String;
    var basisName:String;

    public def this(inpFile:String) {
        this.inputFileName = inpFile;
        try {
            val inp = new JobInput();
            inp.make(inpFile);

            mol = inp.getMolecule();
            basisName = inp.getBasisName();
        } catch(e:Exception) {
            throw new Exception("Unable to read input file: "+inpFile, e);
        }
    } 

    public def runHF() {
        val jd = JobDefaults.getInstance();
        Console.OUT.println("PumjaRasaayani shunya.tri, Quantum Chemistry program in x10, v0.4");
        Console.OUT.println("" + Place.numPlaces() + " places, " + Runtime.NTHREADS + " threads per place");

        mol.transformToSNO(jd.roZ); // Other transformations exist 

        jd.rad=mol.getRadius(jd.roZ);
        Console.OUT.printf("rad/PI=%f\n",jd.rad/3.1415926535);
        if (jd.roOn!=0n && jd.rad>3.1415926535 && jd.maxIterations > 0n) Console.OUT.printf("WARNING: Full-Coulomb RO is not valid for rad>PI\n"); // See line 87 HartreeFockSCFMethod.x10
        if (jd.roOn==0n && jd.roZ!=1.0) Console.OUT.printf("WARNING: Coulomb RO is off: roZ=%f is not necessary.\n",jd.roZ);

        printInput();

        if (jd.useMta) { runMTA(); return; }

        Console.OUT.println("Number of atoms: " + mol.getNumberOfAtoms());
        Console.OUT.println("Number of electrons: " + mol.getNumberOfElectrons());
        val timer = new Timer(3);
        timer.start(0);

        timer.start(1);
        val bsf = new BasisFunctions(mol, basisName, getBasisDirName(inputFileName));
        timer.stop(1);
        Console.OUT.println("Number of basis functions: " + bsf.getBasisFunctions().size());
        Console.OUT.printf("\tTime for setting up basis functions: %.3f seconds\n", (timer.total(1) as Double) / 1e9);
        
        timer.start(2);
        Console.OUT.printf("Computing 1e integrals...\n");
        val oneE = new OneElectronIntegrals(bsf, mol, inputFileName);
        timer.stop(2);
        Console.OUT.printf("\tTime for computing 1E integrals: %.3f seconds\n", (timer.total(2) as Double) / 1e9);
        val hfscf = new HartreeFockSCFMethod(mol, oneE, bsf);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.printf("\n\nTotal time since start: %.3f seconds\n\n", (timer.total(0) as Double) / 1e9);
    }

    private def runHF(fragment:Fragment) {
        val timer = new Timer(3);
        timer.start(0);
        
        Console.OUT.println("\nfragment:");
        Console.OUT.println(fragment);

        timer.start(1);
        val jd = JobDefaults.getInstance();
        val bsf = new BasisFunctions(fragment, basisName, getBasisDirName(inputFileName));
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        timer.stop(1);
        Console.OUT.printf("\tTime for setting up basis functions: %.3f seconds\n\n", (timer.total(1) as Double) / 1e9);

        timer.start(2);
        val oneE = new OneElectronIntegrals(bsf, fragment, inputFileName);
        Console.OUT.println("\nComputed one-electron integrals.");
        timer.stop(2);
        Console.OUT.printf("\tTime for computing 1E integrals: %.3f seconds\n\n", (timer.total(2) as Double) / 1e9);

        val hfscf = new HartreeFockSCFMethod(fragment, oneE, bsf);
        hfscf.scf();
        timer.stop(0);

        Console.OUT.printf("\nTotal time since start: %.3f seconds\n", (timer.total(0) as Double) / 1e9);

        fragment.energy = hfscf.getEnergy();
    }

    public def runMTA() {
        val timer = new Timer(1);
        timer.start(0);

        val fragmentor = new Fragmentor(5.67, 30n);  // TODO, parameters to be taken from user
        
        // first generate the fragments, along with cardinality expression
        val fragments = fragmentor.fragment(mol);

        // run hf for all all the fragments, 
        // TODO: how to parallelize?
        for(fragment in fragments) {
            runHF(fragment); 
        } // end for

        // collect and patch the results using cardinality expression
        var ene:Double = 0.0;
        for(fragment in fragments) {
            ene += fragment.energy * fragment.cardinalitySign;
        } // end for 
        timer.stop(0);

        Console.OUT.printf("Final MTA energy: %.6f a.u.\n", ene);
        Console.OUT.printf("\n-End of MTA run-\n\nTotal time since start: %.3f seconds\n", (timer.total(0) as Double) / 1e9);
    }

    private def printInput() {
        Console.OUT.printf("----------------------------------------------------------\n");
        Console.OUT.println("Input (including defaults):");
        Console.OUT.printf("----------------------------------------------------------\n");

        Console.OUT.println(mol.getName());
        Console.OUT.println("basis " + basisName);
        Console.OUT.println("molecule");
        Console.OUT.println(mol + "end");

        val jd = JobDefaults.getInstance();

        Console.OUT.println("scf_method " + jd.scfMethod);
        Console.OUT.println("scf_max " + jd.maxIterations);
        Console.OUT.println("scf_converge " + jd.energyTolerance);
        Console.OUT.println("diis_start " + jd.diisStartThreshold);
        Console.OUT.println("diis_converge " + jd.diisConvergenceThreshold);
        Console.OUT.println("diis_subspace_size " + jd.diisSubspaceSize);
        Console.OUT.println("guess " + jd.guess);
        if (jd.gMatrixParallelScheme != GMatrix.DEFAULT_GMATTYPE) {
            Console.OUT.println("gmat_parallel_scheme " + jd.gMatrixParallelScheme);
        }
        if (jd.useMta) {
            Console.OUT.println("fragment mta");
        }

        Console.OUT.printf("----------------------------------------------------------\n");
    }

    private static def getBasisDirName(inpFile : String) {
        val inputFile = new File(inpFile);
        val parent = inputFile.getParentFile();
        val inputFileDir = parent == null ? "" : (parent.getPath() + File.SEPARATOR);
        return inputFileDir + "basis";
    }

    public static def main(args:Rail[String]) {
        if (args.size != 1L) {
            Console.ERR.println("usage: pumjarasaayani <inputFile>");
            System.setExitCode(1n);
        } else {
            val qmApp = new PumjaRasaayani(args(0));
            qmApp.runHF();
        }
    }
}

