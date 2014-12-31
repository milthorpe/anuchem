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
import x10x.vector.Point3d;

import au.edu.anu.chem.mm.AtomType;

/**
 * This class reads XYZ molecular structure files of the following format:
 *  <number of atoms>
 *  <title>
 *  species x y z // repeated, positions in nm
 */
public class XYZStructureFileReader { 
    var fileName:String;

    public def this(fileName:String) { 
        this.fileName = fileName;
    }

    public def readParticleData(particleData:AnummParticleData, forceField:ForceField) {
        val atomTypes = forceField.getAtomTypes();
        val file = new FileReader(new File(fileName));
        // line 1: number of atoms
        val numAtoms = Int.parseInt(file.readLine().trim());

        // line 2: structure title
        val title = file.readLine();
        particleData.description = title;

        // lines 3..*: atom positions
        for (i in 0..(numAtoms-1)) {
            val words = file.readLine().split(" ");
            val symbol = words(0);
            var species:Int = -1n;
            for (j in 0..(atomTypes.size-1)) {
                val atomType = atomTypes(j);
                if (symbol.equals(atomType.name)) {
                    species = j as Int;
                    Console.OUT.println("found species " + species + " " + atomType.name + " " + atomType.mass);
                    break;
                }
            }
            if (species == -1n) {
                throw new IllegalArgumentException("no species found for symbol " + symbol);
            }
            val center = Point3d(Double.parseDouble(words(1).trim()),
                                 Double.parseDouble(words(2).trim()),
                                 Double.parseDouble(words(3).trim())) * 10.0;
            particleData.addAtom(i+1, species, center);
        }
        file.close();

    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

