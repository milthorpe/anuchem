/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
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
 *  <title>
 *  <comments>
 *  <comments>
 *  numAtoms [3] ...
 *  x[10] y[10] z[10] species // repeated for number of atoms
 *  <bondType>  <atomIndex1> <atomIndex2> // repeated for number of bonds
 * M  END
 * @see http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
 */
public class MOLStructureFileReader { 
    var fileName : String;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readMolecule(speciesSpecs:Rail[SpeciesSpec]) : Molecule[MMAtom] {
        val file = new FileReader(new File(fileName));
        val title = file.readLine();
        file.readLine(); // line 2 contains file meta-info
        file.readLine(); // line 3 contains comment
        val numAtomsLine = file.readLine();
        val numAtoms = Int.parseInt(numAtomsLine.substring(0,3));
        var molecule : Molecule[MMAtom] = new Molecule[MMAtom](title);
        // following lines up to numAtoms contain atom data
        for(var i:Int=0; i<numAtoms; i++) {
            val line = file.readLine();
            val x = Double.parseDouble(line.substring(0,10));
            val y = Double.parseDouble(line.substring(10,20));
            val z = Double.parseDouble(line.substring(20,30));
            val symbol = line.substring(31,34).trim();
            var species:Int = -1;
            var speciesSpec:SpeciesSpec = null;
            for (j in 0..(speciesSpecs.size-1)) {
                if (symbol.equals(speciesSpecs(j).name)) {
                    species = j as Int;
                    break;
                }
            }
            if (species == -1) {
                throw new IllegalArgumentException("no species found for symbol " + symbol);
            }
            molecule.addAtom(new MMAtom(species, Point3d(x, y, z), speciesSpec.mass, speciesSpec.charge));
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

    public def setFileName(fileName : String) {
        this.fileName = fileName;
    }
}

