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
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.io.File;
import x10.io.FileReader;
import x10.io.IOException;
import au.edu.anu.chem.BondType;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import x10x.vector.Point3d;

/**
 * This class reads Symyx molfiles with the following format:
 *  <number of atoms>
 *  <title>
 *  species x y z // repeated for number of atoms
 *  <bondType>  <atomIndex1> <atomIndex2>
 * @see http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
 */
public class MOLStructureFileReader { 
    var fileName : String;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readMolecule() : Molecule[MMAtom] throws IOException {
        val file = new FileReader(new File(fileName));
        val title = file.readLine();
        file.readLine(); // line 2 contains file meta-info
        file.readLine(); // line 3 contains comment
        val numAtoms = Int.parseInt(file.readLine()); // line 4 contains numAtoms etc.
        var molecule : Molecule[MMAtom]{self.at(this)} = new Molecule[MMAtom](title);
        // following lines up to numAtoms contain atom data
        for(var i:Int=0; i<numAtoms; i++) {
            val line = file.readLine();
            val x = Double.parseDouble(line.substring(0,10));
            val y = Double.parseDouble(line.substring(10,20));
            val z = Double.parseDouble(line.substring(20,30));
            val symbol = line.substring(31,34).trim();
            molecule.addAtom(new MMAtom(symbol, Point3d(x, y, z), 0.0));
        }
        // following lines until "M  END" contain bonds
        var line : String = file.readLine();
        while (line.indexOf("END") < 0) {
            val firstAtomIndex = Int.parseInt(line.substring(0,3))-1;
            val secondAtomIndex = Int.parseInt(line.substring(3,6))-1;
            val bondVal = Int.parseInt(line.substring(6,9));
            var bondType : BondType;
            switch (bondVal) {
                case (1):
                    bondType = BondType.SINGLE_BOND;
                    break;
                case (2):
                    bondType = BondType.DOUBLE_BOND;
                    break;
                case (3):
                    bondType = BondType.TRIPLE_BOND;
                    break;
                case (4):
                    bondType = BondType.AROMATIC_BOND;
                    break;
                default:
                    Console.ERR.println("Bond type " + bondVal + " not found. Defaulting to single bond.");
                    bondType = BondType.SINGLE_BOND;
            }
            val firstAtom = molecule.getAtom(firstAtomIndex);
            val secondAtom = molecule.getAtom(secondAtomIndex);
            //Console.OUT.println("bonding " + firstAtom + " to " + secondAtom + " via " + bondType);
            firstAtom.addBond(bondType, secondAtom);

            line = file.readLine();
        }
        file.close();
        return molecule;
    }

    public safe def setFileName(fileName : String!) {
        this.fileName = fileName;
    }
}

