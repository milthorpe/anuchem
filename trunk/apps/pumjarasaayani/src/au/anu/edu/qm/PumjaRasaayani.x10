/**
 * PumjaRasaayani .. first x10 Quantum Chemistry code (??!!) ;-)
 * Pumja : Sanskrit for Quantum
 * Rasaayn : Chemical (Rasayan Shastra : Chemistry / Chemical Science)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.io.Console;

public class PumjaRasaayani { 
    public def this() {     
    } 

    public def runIt() : void {
        Console.OUT.println("PumjaRasaayani shunya.eak, Quantum Chemisty program in x10, v0.1");

        val mol = new Molecule("h2");

        // H2, 1 a.u. apart
        mol.addAtom(new Atom("H", 0.0, 0.0, 0.0));
        mol.addAtom(new Atom("H", 1.0, 0.0, 0.0));

        Console.OUT.println("\nInput deck:");
        Console.OUT.println(mol);

        val bsf  = new BasisFunctions(mol, "sto3g");
        val oneE = new OneElectronIntegrals(bsf, mol);
        Console.OUT.println("HCore");
        Console.OUT.println(oneE.getHCore());   
        Console.OUT.println("Overlap");
        Console.OUT.println(oneE.getOverlap());   

        val twoE = new TwoElectronIntegrals(bsf, true);

        val hfscf = new HartreeFockSCFMethod(mol, oneE, twoE);
        hfscf.scf();
    }

    public static def main(args:Rail[String]) : void {
        (new PumjaRasaayani()).runIt();
    }
}

