/**
 * AtomInfo.x10
 *
 * Basic information on Atom types.
 *
 * @author: V.Ganesh
 */
package au.edu.anu.chem;

import x10.io.*;
import x10.util.*;

public class AtomInfo {
    var atomicNumbers:HashMap[String, Int]{self.at(this)};
    val atomInfoFile = "AtomInfo.conf";

    var madeIt:Boolean;

    private def this() {
        madeIt = false;
    }

    private def make() {
        if (madeIt) return;

        atomicNumbers = new HashMap[String, Int]();

        /** 
        try {
            val fil = new FileReader(new File(atomInfoFile));
            var lin:String;

            while(true) {
                lin = ""; 
                try {
                   lin = fil.readLine();
                } catch(ignored:EOFException) {  break; }

                x10.io.Console.OUT.println(lin.split(" "));                
            } // end while

            fil.close();
        } catch(e:Exception) {
            x10.io.Console.OUT.println("Error : " + e);
        } // end of try .. catch block
        **/

        atomicNumbers.put("H", 1);
        atomicNumbers.put("He", 2);
        atomicNumbers.put("Li", 3);
        atomicNumbers.put("Be", 4);
        atomicNumbers.put("B", 5);
        atomicNumbers.put("C", 6);
        atomicNumbers.put("N", 7);
        atomicNumbers.put("O", 8);
        atomicNumbers.put("F", 9);
        atomicNumbers.put("Ne", 10);

        madeIt = true;
    }

    private static _theInstance = new AtomInfo();

    public static def getInstance() {
        _theInstance.make();

        return _theInstance;
    }

    public def getAtomicNumber(atm:Atom{self.at(this)}) = atomicNumbers.getOrElse(atm.symbol, -1);
}

