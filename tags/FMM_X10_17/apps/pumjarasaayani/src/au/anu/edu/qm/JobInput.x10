package au.anu.edu.qm;

import x10.io.*;
import x10.util.*;

/**
 * This class expects the input file to be in follwing format
 *
 *  <number of atoms>
 *  <title> <basis>
 *  symbol x y z
 *  ...
 */
public class JobInput { 
    var molecule:Molecule;
    var basisName:String;

    public def this() { } 

    public def make(inpFile:String) throws IOException { 
       readInp(inpFile);
    } 

    private def readInp(inpFile:String) throws IOException { 
       val fil = new FileReader(new File(inpFile));

       val noOfAtoms = Int.parseInt(fil.readLine());
       var words:ArrayList[String] = Utility.split(fil.readLine(), ' ');
       var line:String;
 
       molecule = new Molecule(words(0));
       basisName = words(1);

       for(var i:Int=0; i<noOfAtoms; i++) { 
         line = fil.readLine();         
         words = Utility.split(line, ' ');      

         // for(var j:Int=0; j<words.size(); j++) x10.io.Console.OUT.println(words.get(j));

         molecule.addAtom(new Atom(words.get(0), 
                                   Double.parseDouble(words.get(1)),
                                   Double.parseDouble(words.get(2)),
                                   Double.parseDouble(words.get(3))));
       } // end while              

       fil.close();
    }

    public def getMolecule()  : Molecule = molecule;
    public def getBasisName() : String   = basisName;
}

