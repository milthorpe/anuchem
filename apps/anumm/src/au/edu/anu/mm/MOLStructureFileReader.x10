/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Josh Milthorpe 2010.
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.io.File;
import x10.io.FileReader;
import x10.io.IOException;

import au.edu.anu.chem.BondType;
import au.edu.anu.mm.uff.UniversalForceField;
import x10x.vector.Point3d;

/**
 * This class reads Symyx molfiles with the following format:
 *  <title>
 *  <comments>
 *  <comments>
 *  numAtoms [3] ...
 *  x[10] y[10] z[10] species // repeated for number of atoms, in Angstroms
 *  <atomIndex1> <atomIndex2> <bondType> // repeated for number of bonds
 * M  END
 * @see http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
 */
public class MOLStructureFileReader { 
    var fileName:String;

    public def this(fileName:String) { 
        this.fileName = fileName;
    }

    public def readParticleData(particleData:AnummParticleData, forceField:ForceField) {
        val atomTypes = forceField.getAtomTypes();
        val file = new FileReader(new File(fileName));
        val title = file.readLine();
        file.readLine(); // line 2 contains file meta-info
        file.readLine(); // line 3 contains comment
        val metaLine = file.readLine();
        val numAtoms = Int.parseInt(metaLine.substring(0n,3n).trim());
        val numBonds = Int.parseInt(metaLine.substring(3n,6n).trim());

        particleData.description = title;

        // following lines up to numAtoms contain atom data
        for (i in 0..(numAtoms-1)) {
            val line = file.readLine();
            val center = Point3d(Double.parseDouble(line.substring( 0n,10n).trim()),
                                 Double.parseDouble(line.substring(10n,20n).trim()),
                                 Double.parseDouble(line.substring(20n,30n).trim())) * 0.1; // convert to nm
            val symbol = line.substring(31n,34n).trim();
            var species:Int = -1n;
            for (j in 0..(atomTypes.size-1)) {
                val atomType = atomTypes(j);
                if (symbol.equals(atomType.name)) {
                    species = j as Int;
                    //Console.OUT.println("found species " + species + " " + speciesSpec.name + " " + speciesSpec.mass);
                    break;
                }
            }
            if (species == -1n) {
                throw new IllegalArgumentException("no species found for symbol " + symbol);
            }
            //Console.OUT.println("position " + center);
            particleData.addAtom(i+1, species, center);
        }

        val bonds = new Rail[Bond](numBonds);
        // following lines until "M  END" contain bonds
        var line : String = file.readLine();
        for (i in 0..(numBonds-1)) {
            val atom1Index = Int.parseInt(line.substring(0n,3n).trim())-1;
            val atom2Index = Int.parseInt(line.substring(3n,6n).trim())-1;
            val bondType = Int.parseInt(line.substring(6n,9n).trim());

            bonds(i) = new Bond(atom1Index, atom2Index, forceField.getBondTypeIndex(particleData.atomTypeIndex(atom1Index), particleData.atomTypeIndex(atom2Index), bondType));
            //Console.OUT.println("bond " + atom1Index + " to " + atom2Index + " type " + bondType + " bondTypeIndex " + bonds(i).typeIndex);

            line = file.readLine();
        }
        particleData.bonds = bonds;
        val dummy = new MoleculeType();
        dummy.name = "unknown";
        particleData.moleculeTypes = [ dummy as MoleculeType ];
        particleData.moleculeTypeIndex = [ 0n as Int ];

        particleData.atomTypes = forceField.getAtomTypes();
        particleData.boxEdges = particleData.getMaxExtent();
    
        file.close();
    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

