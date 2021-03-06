/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */
package au.edu.anu.chem.mm;

import x10.io.File;
import x10.io.FileReader;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class reads GROMACS coordinate files of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 *  <box x> <box y> <box z>
 * N.B. GROMACS coordinates are in nm; these are converted to Angstroms for use in ANU-Chem
 * TODO bonding, velocity
 * TODO charges are hardcoded to -0.82 for OW, +0.41 for HW
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsStructureFileReader { 
    var fileName : String;
    public var boxEdges:Vector3d;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readMolecule() : Molecule[MMAtom] {
        val file = new FileReader(new File(fileName));
        val title = file.readLine().split(" ");
        val numAtoms = Int.parseInt(file.readLine());
        var molecule : Molecule[MMAtom] = new Molecule[MMAtom](title(0));
        for(var i:Int=0n; i<numAtoms; i++) {
            val line = file.readLine();
            var species:Int=-1n;
            val atomType = line.substring(10n,15n).trim();
            val charge : Double;
            val mass : Double;
            if (atomType.startsWith("OW")) {
                species = 0n;
                charge = -0.82;
                mass = 15.99491461956;
            } else if (atomType.startsWith("HW")) {
                species = 1n;
                charge = 0.41;
                mass = 1.00794;
            } else {
                Console.ERR.println("Unknown atom type line " + (i+2) + ": " + atomType);
                mass = 0.0;
                charge = 0.0;
            }
            // multiply coords by 10 to convert nm to Angstroms
            molecule.addAtom(new MMAtom(species,
                                         Point3d(Double.parseDouble(line.substring(20n,28n)) * 10.0,
                                                 Double.parseDouble(line.substring(28n,36n)) * 10.0,
                                                 Double.parseDouble(line.substring(36n,44n)) * 10.0
                                         ),
                                        mass,
                                        charge
                        ));
        }

        // read box edge lengths
        val boxLine = file.readLine();
        val x = Double.parseDouble(boxLine.substring( 0n,10n).trim()) * 10.0;
        val y = Double.parseDouble(boxLine.substring(10n,20n).trim()) * 10.0;
        val z = Double.parseDouble(boxLine.substring(20n,30n).trim()) * 10.0;
        this.boxEdges = Vector3d(x,y,z);
        
        file.close();
        return molecule;
    }

    public def setFileName(fileName : String) {
        this.fileName = fileName;
    }
}

