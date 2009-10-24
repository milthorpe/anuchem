/**
 * PumjaRasaayani .. first, fully working, x10 Quantum Chemistry code (??!!) ;-)
 * 
 * Pumja : Sanskrit for Quantum
 * Rasaayn : Chemical (Rasayan Shastra : Chemistry / Chemical Science)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.io.Console;

public class PumjaRasaayani { 
    var mol:Molecule;
    var basisName:String;

    public def this() { }

    public def make() {
        initDefault();
    } 

    public def make(inpFile:String) {     
        try {
          val inp = new JobInput();
          inp.make(inpFile);
        
          mol = inp.getMolecule();
          basisName = inp.getBasisName();
        } catch(e:Exception) {
          Console.OUT.println("Unable to read input file: " + inpFile); 
          Console.OUT.println("Using defaults!");
          initDefault();
        } // end of try .. catch block
    } 

    private def initDefault() { 
        mol = new Molecule("h2");

        // H2, 1 a.u. apart
        mol.addAtom(new Atom("H", 0.0, 0.0, 0.0));
        mol.addAtom(new Atom("H", 1.0, 0.0, 0.0));

        // H2O - differed till HashMap is fixed
        // mol.addAtom(new Atom("O", -0.015283, 0.044743, 6.043909));
        // mol.addAtom(new Atom("H", -0.195295, 0.167197, 6.999322));
        // mol.addAtom(new Atom("H", 0.962423, -0.010198, 6.013442));

        // default basis is sto3g
        basisName = "sto3g";
    }

    public def runIt() {
        Console.OUT.println("PumjaRasaayani shunya.eak, Quantum Chemisty program in x10, v0.1");

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        val bsf  = new BasisFunctions();
        bsf.make(mol, basisName, "basis");
        Console.OUT.println("\nSetting up basis functions over.");
        
        val oneE = new OneElectronIntegrals(bsf, mol);
        oneE.make();
        Console.OUT.println("\nComputed one-electron integrals.");
        Console.OUT.println("HCore");
        Console.OUT.println(oneE.getHCore());   
        Console.OUT.println("Overlap");
        Console.OUT.println(oneE.getOverlap());   

        val twoE = new TwoElectronIntegrals(); 
        twoE.make(bsf, true);
        Console.OUT.println("\nNumber of 2E integrals: " + twoE.getNumberOfIntegrals());
        Console.OUT.println("\nComputed two-electron integrals. If direct, this is skipped for now.");

        val hfscf = new HartreeFockSCFMethod(mol, oneE, twoE);
        hfscf.scf();
    }

    public static def main(args:Rail[String]) {
        val qmApp = new PumjaRasaayani();

        if (args.length == 0) qmApp.make();
        else                  qmApp.make(args(0));
 
        qmApp.runIt();
    }
}
