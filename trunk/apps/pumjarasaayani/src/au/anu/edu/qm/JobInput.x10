package au.anu.edu.qm;

import x10.io.*;
import au.edu.anu.chem.Molecule;
import x10x.vector.Point3d;

/**
 * This class expects the input file to be in follwing format
 *
 *  <number of atoms>
 *  <title> <basis>
 *  symbol x y z
 *  ...
 */
public class JobInput { 
    global var molecule:Molecule[QMAtom]{self.at(this)};
    var basisName:String;

    public def this() { } 

    public def make(inpFile:String) throws IOException { 
       readInp(inpFile);
    } 

    private def readInp(inpFile:String) throws IOException { 
       val fil = new FileReader(new File(inpFile));

       val noOfAtoms = Int.parseInt(fil.readLine());
       val words = fil.readLine().split(" ");
 
       molecule = new Molecule[QMAtom]();
       molecule.make(words(0));
       basisName = words(1);

       for(var i:Int=0; i<noOfAtoms; i++) { 
         val wrd = fil.readLine().split(" ");         

         // for(var j:Int=0; j<words.size(); j++) x10.io.Console.OUT.println(words.get(j));

         molecule.addAtom(new QMAtom(wrd(0), 
                                   new Point3d(Double.parseDouble(wrd(1)),
                                               Double.parseDouble(wrd(2)),
                                               Double.parseDouble(wrd(3))
                                    )
                            ));
       } // end while              

       fil.close();
    }

    public def getMolecule()  : Molecule[QMAtom]{self.at(this)} = molecule;
    public def getBasisName() : String   = basisName;
}

