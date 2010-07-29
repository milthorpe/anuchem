/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import x10.io.*;
import au.edu.anu.chem.Molecule;
import x10x.vector.Point3d;

/**
 * This class expects the input file to be in following format.
 * Coordinates are expected to be in a.u.
 *
 *  <number of atoms>
 *  <title> <basis>
 *  symbol x y z
 *  ...
 */
public class JobInput { 
    var molecule:Molecule[QMAtom]{self.at(this)};
    var basisName:String;

    public def this() { } 

    public def make(inpFile:String) throws IOException { 
       readInp(inpFile);
    } 

    private def readInp(inpFile:String) throws IOException { 
       x10.io.Console.OUT.println("NOTE: the input file is expected to be in ** a.u. ** units");

       val fil = new FileReader(new File(inpFile));

       val noOfAtoms = Int.parseInt(fil.readLine());
       val words = fil.readLine().split(" ");
 
       molecule = new Molecule[QMAtom](words(0));
       basisName = words(1);

       for(var i:Int=0; i<noOfAtoms; i++) { 
         val wrd = fil.readLine().split(" ");         

         molecule.addAtom(new QMAtom(wrd(0), 
                                       Point3d(Double.parseDouble(wrd(1)),
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

