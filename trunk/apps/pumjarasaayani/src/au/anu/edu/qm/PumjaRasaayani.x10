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
import au.edu.anu.util.Timer;
import x10x.vector.Point3d;

import au.anu.edu.qm.mta.Fragment; 
import au.anu.edu.qm.mta.Fragmentor; 

public class PumjaRasaayani { 
    var mol:Molecule[QMAtom]{self.at(this)};
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
        if (isMTA) { runMTA(); return; } // is it an MTA run?
 
        val timer = new Timer(3);
        timer.start(0);

        Console.OUT.println("PumjaRasaayani shunya.dau, Quantum Chemisty program in x10, v0.2");

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);

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

        val twoE = new TwoElectronIntegrals(bsf, mol, true);
        Console.OUT.println("\nNumber of 2E integrals: " + twoE.getNumberOfIntegrals());
        Console.OUT.println("\nComputed two-electron integrals. If direct, this is skipped for now.");
        Console.OUT.println("Is Direct: " + twoE.isDirect());

        val hfscf = new HartreeFockSCFMethod(mol, oneE, twoE, gMatType);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.println ("\n-End of SCF-\n\nTotal time since start: " + (timer.total(0) as Double) / 1e9 + " seconds\n");
        
        energy = hfscf.getEnergy();
        time   = (timer.total(0) as Double) / 1e9;
    }

    private def runHF(fragment:Fragment!) {
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

        val twoE = new TwoElectronIntegrals(bsf, fragment, true);
        Console.OUT.println("\nNumber of 2E integrals: " + twoE.getNumberOfIntegrals());
        Console.OUT.println("\nComputed two-electron integrals. If direct, this is skipped for now.");
        Console.OUT.println("Is Direct: " + twoE.isDirect());

        val hfscf = new HartreeFockSCFMethod(fragment, oneE, twoE, gMatType);
        hfscf.scf();
        timer.stop(0);
        Console.OUT.println ("\n-End of SCF-\n\nTotal time since start: " + (timer.total(0) as Double) / 1e9 + " seconds\n");

        fragment.energy(hfscf.getEnergy());
        time   = (timer.total(0) as Double) / 1e9;
    }

    public def runMTA() {
        val fragmentor = new Fragmentor(5.67, 30);  // TODO, parameters to be taken from user
        
        // first generate the main fragments 
        val mainFragments = fragmentor.fragment(mol);

        // then generate the cardinality expression

        // run hf for all all the fragments
        for(fragment in mainFragments) {
            runHF(fragment as Fragment!); 
        } // end for

        // collect and patch the results using cardinality expression
    }

    public static def main(args:Rail[String]!) {
        val qmApp = args.length == 0 ? new PumjaRasaayani() : 
                    args.length == 1 ? new PumjaRasaayani(args(0)) : 
                    args.length == 2 ? new PumjaRasaayani(args(0), Int.parseInt(args(1))) : 
                                       new PumjaRasaayani(args(0), Int.parseInt(args(1)), args(2));
        qmApp.runIt();
    }
}

