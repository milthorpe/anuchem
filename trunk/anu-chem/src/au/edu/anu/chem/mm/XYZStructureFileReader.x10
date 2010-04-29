package au.edu.anu.chem.mm;

import x10.io.*;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import x10x.vector.Point3d;

/**
 * This class reads XYZ molecular structure files of the following format:
 *  <number of atoms>
 *  <title>
 *  species x y z // repeated
 */
public class XYZStructureFileReader { 
    var fileName : String;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readMolecule() : Molecule[MMAtom] throws IOException {
        val file = new FileReader(new File(fileName));
        val numAtoms = Int.parseInt(file.readLine());
        val title = file.readLine().split(" ");
        var molecule : Molecule[MMAtom]{self.at(this)} = new Molecule[MMAtom](title(0));
        for(var i:Int=0; i<numAtoms; i++) {
            val words = file.readLine().split(" ");
            molecule.addAtom(new MMAtom(words(0),
                                         Point3d(Double.parseDouble(words(1)),
                                                 Double.parseDouble(words(2)),
                                                 Double.parseDouble(words(3))
                                         ),
                                        0.0,
                                        0.0
                        ));
        }
       file.close();
       return molecule;
    }

    public safe def setFileName(fileName : String!) {
        this.fileName = fileName;
    }
}

