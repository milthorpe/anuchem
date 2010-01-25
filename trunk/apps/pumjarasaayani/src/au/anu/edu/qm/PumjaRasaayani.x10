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
    global var mol:Molecule{self.at(this)};
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
        mol = new Molecule();
        mol.make("h2");

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

        val bsf:BasisFunctions{self.at(this)}  = new BasisFunctions();
        bsf.make(mol, basisName, "basis");
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        
        val oneE:OneElectronIntegrals{self.at(this)} = new OneElectronIntegrals();
        oneE.make(bsf, mol);
        Console.OUT.println("\nComputed one-electron integrals.");
        Console.OUT.println("HCore");
        Console.OUT.println(oneE.getHCore());   
        Console.OUT.println("Overlap");
        Console.OUT.println(oneE.getOverlap());   

        val twoE:TwoElectronIntegrals{self.at(this)} = new TwoElectronIntegrals(); 
        twoE.make(bsf, mol, true);
        Console.OUT.println("\nNumber of 2E integrals: " + twoE.getNumberOfIntegrals());
        Console.OUT.println("\nComputed two-electron integrals. If direct, this is skipped for now.");
        Console.OUT.println("Is Direct: " + twoE.isDirect());

        val hfscf:HartreeFockSCFMethod{self.at(this)} = new HartreeFockSCFMethod();
        hfscf.make(mol, oneE, twoE);
        hfscf.scf();
    }

    public def make(ag:Rail[String]) {
        val args:Rail[String]{self.at(this)} = ag as Rail[String]{self.at(this)}; 

        if (args.length == 0) make();
        else                  make(args(0));
    }

    public static def main(args:Rail[String]) {
        val qmApp   = new PumjaRasaayani();

        qmApp.make(args);
        qmApp.runIt();
    }
}

