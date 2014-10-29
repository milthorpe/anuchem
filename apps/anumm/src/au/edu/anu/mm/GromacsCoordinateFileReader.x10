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
import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class reads GROMACS coordinate files (*.gro) of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 *  <box x> <box y> <box z>
 * TODO bonding, velocity
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsCoordinateFileReader { 
    var fileName : String;
    public var boxEdges:Vector3d;

    public def this(fileName : String) { 
        this.fileName = fileName;
    }

    public def readParticleData(particleData:ParticleData, forceField:ForceField) {
        val speciesSpecs = forceField.getSpeciesSpecs();
        val file = new FileReader(new File(fileName));
        val title = file.readLine(); // TODO should read from topology file
        val numAtoms = Int.parseInt(file.readLine());
        
        particleData.allocateAtoms(numAtoms);
        particleData.description = title;

        for(var i:Int=0n; i<numAtoms; i++) {
            val line = file.readLine();

            val atomType = line.substring(10n,15n).trim();
            // TODO look up atom type from atom name in residues
            var species:Int = -1n;
            for (j in 0..(speciesSpecs.size-1)) {
                val speciesSpec = speciesSpecs(j);
                if (atomType.equals(speciesSpec.name)) {
                    species = j as Int;
                    //Console.OUT.println("found species " + species + " " + speciesSpec.name + " " + speciesSpec.mass);
                    break;
                }
            }
            if (species == -1n) {
                throw new IllegalArgumentException("no species found for symbol " + atomType);
            }
            particleData.species(i) = species;

            particleData.globalIndex(i) = Long.parseLong(line.substring(15n,20n).trim());
            // multiply coords by 10 to convert nm to Angstroms
            particleData.x(i) = Point3d(Double.parseDouble(line.substring(20n,28n)) * 10.0,
                                        Double.parseDouble(line.substring(28n,36n)) * 10.0,
                                        Double.parseDouble(line.substring(36n,44n)) * 10.0);

            if (line.length() > 67) {
                // velocities are included
                particleData.dx(i) = Vector3d(Double.parseDouble(line.substring(44n,52n)) * 10.0,
                                              Double.parseDouble(line.substring(52n,60n)) * 10.0,
                                              Double.parseDouble(line.substring(60n,68n)) * 10.0);
            }
        }

        // read box edge lengths
        val boxLine = file.readLine();
        val x = Double.parseDouble(boxLine.substring( 0n,10n).trim()) * 10.0;
        val y = Double.parseDouble(boxLine.substring(10n,20n).trim()) * 10.0;
        val z = Double.parseDouble(boxLine.substring(20n,30n).trim()) * 10.0;
        this.boxEdges = Vector3d(x,y,z);
        
        file.close();
    }

    public def setFileName(fileName : String) {
        this.fileName = fileName;
    }
}

