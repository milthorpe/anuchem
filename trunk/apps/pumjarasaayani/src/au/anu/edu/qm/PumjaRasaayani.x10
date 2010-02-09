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
import au.edu.anu.chem.Molecule;
import x10x.vector.Point3d;

public class PumjaRasaayani { 
    global var mol:Molecule[QMAtom]{self.at(this)};
    var basisName:String;

    public def this() {
        initDefault();
    } 

    public def this(inpFile:String) {     
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
        mol = new Molecule[QMAtom]();
        mol.make("h2");

        // H2, 1 a.u. apart
        mol.addAtom(new QMAtom("H", new Point3d(0.0, 0.0, 0.0)));
        mol.addAtom(new QMAtom("H", new Point3d(1.0, 0.0, 0.0)));

        // default basis is sto3g
        basisName = "sto3g";
    }

    public def runIt() {
        Console.OUT.println("PumjaRasaayani shunya.eak, Quantum Chemisty program in x10, v0.1");

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);

        Console.OUT.println("\nSetting up basis set: " + basisName);

        val bsf = new BasisFunctions(mol, basisName, "basis");
        Console.OUT.println("\nUsing " + bsf.getBasisFunctions().size() + " basis functions.");
        
        val oneE = new OneElectronIntegrals(bsf, mol);
        Console.OUT.println("\nComputed one-electron integrals.");
        // Console.OUT.println("HCore");
        // Console.OUT.println(oneE.getHCore());   
        // Console.OUT.println("Overlap");
        // Console.OUT.println(oneE.getOverlap());   

        val twoE = new TwoElectronIntegrals(bsf, mol, true);
        Console.OUT.println("\nNumber of 2E integrals: " + twoE.getNumberOfIntegrals());
        Console.OUT.println("\nComputed two-electron integrals. If direct, this is skipped for now.");
        Console.OUT.println("Is Direct: " + twoE.isDirect());

        val hfscf = new HartreeFockSCFMethod(mol, oneE, twoE);
        hfscf.scf();
    }

    public static def main(args:Rail[String]!) {
        val qmApp = args.length == 0 ? new PumjaRasaayani() : new PumjaRasaayani(args(0));
        qmApp.runIt();
    }
}

