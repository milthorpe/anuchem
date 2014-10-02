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
 *  <atomIndex1> <atomIndex2> <bondType> // repeated for number of bonds
 * M  END
 * @see http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
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
        val numAtoms = Int.parseInt(numAtomsLine.substring(0n,3n).trim());
        var molecule : Molecule[MMAtom] = new Molecule[MMAtom](title);
        // following lines up to numAtoms contain atom data
        for(var i:Long=0; i<numAtoms; i++) {
            val line = file.readLine();
            val x = Double.parseDouble(line.substring(0n,10n).trim());
            val y = Double.parseDouble(line.substring(10n,20n).trim());
            val z = Double.parseDouble(line.substring(20n,30n).trim());
            val symbol = line.substring(31n,34n).trim();
            var species:Int = -1n;
            var speciesSpec:SpeciesSpec = null;
            for (j in 0..(speciesSpecs.size-1)) {
                speciesSpec = speciesSpecs(j);
                if (speciesSpec != null && symbol.equals(speciesSpec.name)) {
                    species = j as Int;
                    break;
                }
            }
            if (species == -1n) {
                throw new IllegalArgumentException("no species found for symbol " + symbol);
            }
            molecule.addAtom(new MMAtom(species, Point3d(x, y, z), speciesSpec.mass, speciesSpec.charge));
        }
        // following lines until "M  END" contain bonds
        var line : String = file.readLine();
        while (line.indexOf("END") < 0) {
            val firstAtomIndex = Int.parseInt(line.substring(0n,3n).trim())-1;
            val secondAtomIndex = Int.parseInt(line.substring(3n,6n).trim())-1;
            val bondVal = Int.parseInt(line.substring(6n,9n).trim());
            var bondType : BondType;
            switch (bondVal) {
                case (1n):
                    bondType = BondType.SINGLE_BOND;
                    break;
                case (2n):
                    bondType = BondType.DOUBLE_BOND;
                    break;
                case (3n):
                    bondType = BondType.TRIPLE_BOND;
                    break;
                case (4n):
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

