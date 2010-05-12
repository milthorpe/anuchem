/**
 * AtomInfo.x10
 *
 * Basic information on Atom types.
 * Only H to Ne filled at the moment.
 *
 * @author: V.Ganesh
 */
package au.edu.anu.chem;

import x10.io.*;
import x10.util.*;

public class AtomInfo {
    /** atomic numbers */
    var atomicNumbers:HashMap[String, Int]{self.at(this)};
   
    /** covalent radii, in a.u. */
    var covalentRadii:HashMap[String, Double]{self.at(this)};

    /** vdW radii, in a.u. */
    var vdWRadii:HashMap[String, Double]{self.at(this)};

    var madeIt:Boolean;

    private def this() {
        madeIt = false;
    }

    private def make() {
        if (madeIt) return;

        atomicNumbers = new HashMap[String, Int]();

        atomicNumbers.put("H",  1);
        atomicNumbers.put("He", 2);
        atomicNumbers.put("Li", 3);
        atomicNumbers.put("Be", 4);
        atomicNumbers.put("B",  5);
        atomicNumbers.put("C",  6);
        atomicNumbers.put("N",  7);
        atomicNumbers.put("O",  8);
        atomicNumbers.put("F",  9);
        atomicNumbers.put("Ne", 10);

        covalentRadii = new HashMap[String, Double]();
        covalentRadii.put("H",  0.4346);
        covalentRadii.put("He", 0.6047);
        covalentRadii.put("Li", 0.6047);
        covalentRadii.put("Be", 1.7008);
        covalentRadii.put("B",  1.5496);
        covalentRadii.put("C",  1.4551);
        covalentRadii.put("N",  1.4173);
        covalentRadii.put("O",  1.3228);
        covalentRadii.put("F",  1.3417);
        covalentRadii.put("Ne", 1.3417);

        vdWRadii = new HashMap[String, Double]();
        vdWRadii.put("H",  2.2677);
        vdWRadii.put("He", 2.6456);
        vdWRadii.put("Li", 2.6456);
        vdWRadii.put("Be", 2.5936);
        vdWRadii.put("B",  1.5023);
        vdWRadii.put("C",  3.2125);
        vdWRadii.put("N",  2.9291);
        vdWRadii.put("O",  2.8724);
        vdWRadii.put("F",  2.7779);
        vdWRadii.put("Ne", 2.7779);

        madeIt = true;
    }

    private static _theInstance = new AtomInfo();

    public static def getInstance() {
        _theInstance.make();

        return _theInstance;
    }

    public def getAtomicNumber(atm:Atom{self.at(this)}) = atomicNumbers.getOrElse(atm.symbol, -1);
    public def getCovalentRadius(atm:Atom{self.at(this)}) = covalentRadii.getOrElse(atm.symbol, -1);
    public def getVdwRadius(atm:Atom{self.at(this)}) = vdWRadii.getOrElse(atm.symbol, -1);
}

