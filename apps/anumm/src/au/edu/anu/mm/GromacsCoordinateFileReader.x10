/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Josh Milthorpe 2010-2013.
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.io.File;
import x10.io.FileReader;
import x10.util.ArrayList;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class reads GROMACS coordinate files (*.gro) of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 *  <box x> <box y> <box z>
 * N.B. GROMACS coordinates are in nm; these are converted to Angstroms for use in ANU-Chem
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsCoordinateFileReader { 
    var fileName:String;

    public def this(fileName:String) { 
        this.fileName = fileName;
    }

    public def readParticleData(particleData:ParticleData, forceField:ForceField) {
        val file = new FileReader(new File(fileName));
        val title = file.readLine(); // TODO should read from topology file
        val numAtoms = Int.parseInt(file.readLine());
        
        particleData.allocateAtoms(numAtoms);
        particleData.description = title;

        val bonds = new ArrayList[Bond]();

        for(var i:Int=0n; i<numAtoms; i++) {
            val line = file.readLine();
//Console.OUT.println("reading atom: " + line);
            val residueNumber = Int.parseInt(line.substring(0n,5n).trim());
            particleData.residueNumber(i) = residueNumber;

            val residueName = line.substring(5n,10n).trim();

            for (j in 0..(particleData.moleculeTypes.size-1)) {
                if (residueName.equals(particleData.moleculeTypes(j).name)) {
                    particleData.moleculeTypeIndex(residueNumber) = j as Int;
                    break;
                }
            }
            val moleculeType = particleData.moleculeTypes(particleData.moleculeTypeIndex(residueNumber));
//Console.OUT.println("found residue: " + residueName + " residueNumber " + residueNumber + " moleculeType " + moleculeType.name);

            val atomName = line.substring(10n,15n).trim();
            var atomTypeIndex:Int = -1n;
            for (j in 0..(moleculeType.atomTypes.size()-1)) {
                if (atomName.equals(moleculeType.atomNames(j))) {
                    atomTypeIndex = moleculeType.atomTypes(j);
                    //Console.OUT.println("found atomName " + atomName + " atomTypeIndex " + atomTypeIndex);
                    break;
                }
            }
            if (atomTypeIndex == -1n) {
                throw new IllegalArgumentException("no atomTypeIndex found for atomName " + atomName);
            }
            particleData.atomTypeIndex(i) = atomTypeIndex;

            particleData.globalIndex(i) = Long.parseLong(line.substring(15n,20n).trim());
            // multiply coords by 10 to convert nm to Angstroms
            particleData.x(i) = Point3d(Double.parseDouble(line.substring(20n,28n)),
                                        Double.parseDouble(line.substring(28n,36n)),
                                        Double.parseDouble(line.substring(36n,44n)));

            if (line.length() >= 67) {
                // velocities are included
                particleData.dx(i) = Vector3d(Double.parseDouble(line.substring(44n,52n)),
                                              Double.parseDouble(line.substring(52n,60n)),
                                              Double.parseDouble(line.substring(60n,68n)));
            }

            createBondsForParticle(particleData, i, bonds);
        }

        particleData.bonds = bonds.toRail();

        // read box edge lengths
        val boxLine = file.readLine();
        val x = Double.parseDouble(boxLine.substring( 0n,10n).trim());
        val y = Double.parseDouble(boxLine.substring(10n,20n).trim());
        val z = Double.parseDouble(boxLine.substring(20n,30n).trim());
        particleData.boxEdges = Vector3d(x,y,z);
        
        file.close();
    }

    private def createBondsForParticle(particleData:ParticleData, atomNumber:Long, bonds:ArrayList[Bond]) {
        //val moleculeType = particleData.moleculeTypes(particleData.moleculeTypeIndex(atomNumber));
        //val atomTypeIndex = particleData.atomTypeIndex(atomNumber);
        // TODO store bond template for molecule type when reading topology; use here
    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

