/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.io.File;
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
    var gMatType:Int;
    var isMTA:Boolean;

    public def this(inpFile:String, gMatType:Int) {
        this.inputFileName = inpFile;
        this.gMatType = gMatType;
        try {
            val inp = new JobInput();
            inp.make(inpFile);

            mol = inp.getMolecule();
            basisName = inp.getBasisName();
        } catch(e:Exception) {
            throw new Exception("Unable to read input file: "+inpFile, e);
        }
    } 

    public def this(inpFile:String, gMatType:Int, mtaOpt:String) {
        this(inpFile, gMatType);

        this.isMTA = (mtaOpt.equals("-mta"));
    }

    public def runIt() {
        Console.OUT.println("PumjaRasaayani shunya.tri, Quantum Chemistry program in x10, v0.4");

        Console.OUT.println("No. of places: " + Place.MAX_PLACES);
        Console.OUT.println("No. of threads per place: " + Runtime.NTHREADS);

        mol.transformToSNO();

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);
        Console.OUT.println("Number of atoms: " + mol.getNumberOfAtoms());

        if (isMTA) { runMTA(); return; } // is it an MTA run?
 
        val timer = new Timer(3);
        timer.start(0);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        timer.start(1);
        val bsf = new BasisFunctions(mol, basisName, getBasisDirName(inputFileName));
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        timer.stop(1);
        Console.OUT.printf("    Time for setting up basis functions: %.3g milliseconds\n\n", (timer.total(1) as Double) / 1e6);
        
        timer.start(2);
        val oneE = new OneElectronIntegrals(bsf, mol);
        timer.stop(2);
        Console.OUT.printf("    Time for computing 1E integrals: %.3g seconds\n\n", (timer.total(2) as Double) / 1e9);
        // Console.OUT.println("HCore");
        // Console.OUT.println(oneE.getHCore());   
        // Console.OUT.println("Overlap");
        // Console.OUT.println(oneE.getOverlap());   

        val hfscf = new HartreeFockSCFMethod(mol, oneE, bsf, gMatType);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.printf("\n\nTotal time since start: %.3g seconds\n\n", (timer.total(0) as Double) / 1e9);
    }

    private def runHF(fragment:Fragment) {
        val timer = new Timer(3);
        timer.start(0);

        Console.OUT.println("\nFragment Input deck:");
        Console.OUT.println(fragment);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        timer.start(1);
        val bsf = new BasisFunctions(fragment, basisName, getBasisDirName(inputFileName));
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        timer.stop(1);
        Console.OUT.printf("\tTime for setting up basis functions: %.3g milliseconds\n\n", (timer.total(1) as Double) / 1e9);

        timer.start(2);
        val oneE = new OneElectronIntegrals(bsf, fragment);
        Console.OUT.println("\nComputed one-electron integrals.");
        timer.stop(2);
        Console.OUT.printf("\tTime for computing 1E integrals: %.3g seconds\n\n", (timer.total(2) as Double) / 1e9);
        // Console.OUT.println("HCore");
        // Console.OUT.println(oneE.getHCore());
        // Console.OUT.println("Overlap");
        // Console.OUT.println(oneE.getOverlap());

        val hfscf = new HartreeFockSCFMethod(fragment, oneE, bsf, gMatType);
        hfscf.scf();
        timer.stop(0);

        Console.OUT.printf("\n\nTotal time since start: %.3g seconds\n", (timer.total(0) as Double) / 1e9);

        fragment.energy = hfscf.getEnergy();
    }

    public def runMTA() {
        val timer = new Timer(1);
        timer.start(0);

        val fragmentor = new Fragmentor(5.67, 30);  // TODO, parameters to be taken from user
        
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
        Console.OUT.printf("\n-End of MTA run-\n\nTotal time since start: %.3g seconds\n", (timer.total(0) as Double) / 1e9);
    }

    private static def getBasisDirName(inpFile : String) {
        val inputFile = new File(inpFile);
        val parent = inputFile.getParentFile();
        val inputFileDir = parent == null ? "" : (parent.getPath() + File.SEPARATOR);
        return inputFileDir + "basis";
    }

    public static def main(args : Array[String](1)) {
        if (args.size < 1) {
            Console.ERR.println("usage: pumjarasaayani <inputFile> [<gMatType>] [-mta]");
            System.setExitCode(1);
        } else {
            val qmApp = args.size == 1 ? new PumjaRasaayani(args(0), GMatrix.DEFAULT_GMATTYPE) : 
                        args.size == 2 ? new PumjaRasaayani(args(0), Int.parseInt(args(1))) : 
                                         new PumjaRasaayani(args(0), Int.parseInt(args(1)), args(2));
            qmApp.runIt();
        }
    }
}

