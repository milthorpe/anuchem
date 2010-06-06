/*
 * This file is part of ANU Molecular Mechanics (ANUMM).
 * ANUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUMM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

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

