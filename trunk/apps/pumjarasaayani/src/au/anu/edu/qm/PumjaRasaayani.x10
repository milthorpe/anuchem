/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import x10.io.Console;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;
import x10x.vector.Point3d;

import au.anu.edu.qm.mta.Fragment; 
import au.anu.edu.qm.mta.Fragmentor; 

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
    var basisName:String;
    var gMatType:Int;
    var isMTA:Boolean;

    public def this() {
        initDefault();
    } 

    public def this(inpFile:String) {     
        try {
          val inp = new JobInput();
          inp.make(inpFile);
        
          mol = inp.getMolecule();
          basisName = inp.getBasisName();

          gMatType = 0;
        } catch(e:Exception) {
          Console.OUT.println("Unable to read input file: " + inpFile); 
          Console.OUT.println("Using defaults!");
          initDefault();
        } // end of try .. catch block
    }

    public def this(inpFile:String, gMatType:Int) {
        this(inpFile);
        this.gMatType = gMatType;
    } 

    public def this(inpFile:String, gMatType:Int, mtaOpt:String) {
        this(inpFile);
        this.gMatType = gMatType;

        if (mtaOpt.equals("-mta")) this.isMTA = true;
        else                       this.isMTA = false;
    }

    private def initDefault() { 
        mol = new Molecule[QMAtom]("h2");

        // H2, 1 a.u. apart
        mol.addAtom(new QMAtom("H", Point3d(0.0, 0.0, 0.0)));
        mol.addAtom(new QMAtom("H", Point3d(1.0, 0.0, 0.0)));

        // default basis is sto3g
        basisName = "sto3g";
    }

    var energy:Double = 0.0;
    var time:Double = 0.0;

    public def getEnergy() = energy;
    public def getTime() = time;

    public def runIt() {
        Console.OUT.println("PumjaRasaayani shunya.tri, Quantum Chemisty program in x10, v0.3");

        Console.OUT.println("No. of places: " + Place.MAX_PLACES);
        Console.OUT.println("No. of threads per place: " + Runtime.NTHREADS);

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);
        Console.OUT.println("Number of atoms: " + mol.getNumberOfAtoms());

        if (isMTA) { runMTA(); return; } // is it an MTA run?
 
        val timer = new Timer(3);
        timer.start(0);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        timer.start(1);
        val bsf = new BasisFunctions(mol, basisName, "basis");
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        timer.stop(1);
        Console.OUT.println ("\tTime for setting up basis functions: " + (timer.total(1) as Double) / 1e9 + " seconds\n");
        
        timer.start(2);
        val oneE = new OneElectronIntegrals(bsf, mol);
        Console.OUT.println("\nComputed one-electron integrals.");
        timer.stop(2);
        Console.OUT.println ("\tTime for computing 1E integrals: " + (timer.total(2) as Double) / 1e9 + " seconds\n");
        // Console.OUT.println("HCore");
        // Console.OUT.println(oneE.getHCore());   
        // Console.OUT.println("Overlap");
        // Console.OUT.println(oneE.getOverlap());   

        val hfscf = new HartreeFockSCFMethod(mol, oneE, bsf, gMatType);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.println ("\n-End of SCF-\n\nTotal time since start: " + (timer.total(0) as Double) / 1e9 + " seconds\n");
        
        energy = hfscf.getEnergy();
        time   = (timer.total(0) as Double) / 1e9;
    }

    private def runHF(fragment:Fragment) {
        val timer = new Timer(3);
        timer.start(0);

        Console.OUT.println("\nFragment Input deck:");
        Console.OUT.println(fragment);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        timer.start(1);
        val bsf = new BasisFunctions(fragment, basisName, "basis");
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        timer.stop(1);
        Console.OUT.println ("\tTime for setting up basis functions: " + (timer.total(1) as Double) / 1e9 + " seconds\n");

        timer.start(2);
        val oneE = new OneElectronIntegrals(bsf, fragment);
        Console.OUT.println("\nComputed one-electron integrals.");
        timer.stop(2);
        Console.OUT.println ("\tTime for computing 1E integrals: " + (timer.total(2) as Double) / 1e9 + " seconds\n");
        // Console.OUT.println("HCore");
        // Console.OUT.println(oneE.getHCore());
        // Console.OUT.println("Overlap");
        // Console.OUT.println(oneE.getOverlap());

        val hfscf = new HartreeFockSCFMethod(fragment, oneE, bsf, gMatType);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.println ("\n-End of SCF-\n\nTotal time since start: " + (timer.total(0) as Double) / 1e9 + " seconds\n");

        fragment.energy = hfscf.getEnergy();
        time   = (timer.total(0) as Double) / 1e9;
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

        Console.OUT.println("Final MTA energy : " + ene);
        Console.OUT.println ("\n-End of MTA run-\n\nTotal time since start: " + (timer.total(0) as Double) / 1e9 + " seconds\n");
    }

    public static def main(args : Array[String](1)) {
        val qmApp = args.size == 0 ? new PumjaRasaayani() : 
                    args.size == 1 ? new PumjaRasaayani(args(0)) : 
                    args.size == 2 ? new PumjaRasaayani(args(0), Int.parseInt(args(1))) : 
                                     new PumjaRasaayani(args(0), Int.parseInt(args(1)), args(2));
        qmApp.runIt();
    }
}

