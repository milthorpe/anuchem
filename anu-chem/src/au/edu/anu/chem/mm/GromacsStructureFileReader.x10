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
package au.edu.anu.chem.mm;

import x10.io.File;
import x10.io.FileReader;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import x10x.vector.Point3d;

/**
 * This class reads GROMACS coordinate files of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 * N.B. GROMACS coordinates are in nm; these are converted to Angstroms for use in ANU-Chem
 * TODO bonding, velocity
 * TODO charges are hardcoded to -0.82 for OW, +0.41 for HW
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsStructureFileReader { 
    var fileName : String;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readMolecule() : Molecule[MMAtom] {
        val file = new FileReader(new File(fileName));
        val title = file.readLine().split(" ");
        val numAtoms = Int.parseInt(file.readLine());
        var molecule : Molecule[MMAtom] = new Molecule[MMAtom](title(0));
        for(var i:Int=0; i<numAtoms; i++) {
            val line = file.readLine();
            val atomType = line.substring(10,15).trim();
            val charge : Double;
            if (atomType.startsWith("OW")) {
                charge = -0.82;
            } else if (atomType.startsWith("HW")) {
                charge = 0.41;
            } else {
                Console.ERR.println("Unknown atom type line " + (i+2) + ": " + atomType);
                charge = 0.0;
            }
            // multiply coords by 10 to convert nm to Angstroms
            molecule.addAtom(new MMAtom(atomType,
                                         Point3d(Double.parseDouble(line.substring(20,28)) * 10.0,
                                                 Double.parseDouble(line.substring(28,36)) * 10.0,
                                                 Double.parseDouble(line.substring(36,44)) * 10.0
                                         ),
                                        charge
                        ));
        }
       file.close();
       return molecule;
    }

    public def setFileName(fileName : String) {
        this.fileName = fileName;
    }
}

