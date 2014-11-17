/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.io.File;
import x10.io.FileWriter;
import x10.io.Printer;
import x10.util.StringBuilder;

/**
 * This class reads GROMACS coordinate files (*.gro) of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 *  <box x> <box y> <box z>
 * N.B. GROMACS coordinates are in nm; these are converted to Angstroms for use in ANU-Chem
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsCoordinateFileWriter { 
    val fileName:String;
    val out:Printer;

    public def this(fileName:String) { 
        this.fileName = fileName;
        this.out = new File(fileName).printer();
    }

    public def writePositions(particleDataPlh:PlaceLocalHandle[ParticleData], forceField:ForceField) {
        out.println(particleDataPlh().description);
        out.println(particleDataPlh().numAtoms);
        finish for(place in Place.places()) {
            val placeString = at(place) {
                val particleData = particleDataPlh();
                val posString = new StringBuilder();
                for (i in 0..(particleData.numAtoms-1)) {
                    val moleculeTypeIndex = particleData.moleculeTypeIndex(particleData.residueNumber(i));
                    val residueName = particleData.moleculeTypes(moleculeTypeIndex).name;
                    val name = particleData.atomTypes(particleData.atomTypeIndex(i)).name;
                    posString.add(String.format("%5d%3s  %5s%5d", [particleData.residueNumber(i), residueName, name, particleData.globalIndex(i)]));
                    val pos = particleData.x(i);
                    val vel = particleData.dx(i);
                    posString.add(String.format("%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", [pos.i as Any, pos.j, pos.k, vel.i, vel.j, vel.k]));
                }
                posString.toString()
            };
            out.print(placeString);
        }
        val edges = particleDataPlh().boxEdges;
        out.printf("%10.5f%10.5f%10.5f\n", edges.i, edges.j, edges.k);
    }
}

