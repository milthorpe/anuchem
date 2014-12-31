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
import x10.util.HashMap;

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

    public def readParticleData(particleData:AnummParticleData) {
        val file = new FileReader(new File(fileName));
        val title = file.readLine(); // TODO should read from topology file
        val numAtoms = Int.parseInt(file.readLine());
        
        particleData.description = title;

        // maps residue name, number, atom name to globalIndex
        val particleMap = new HashMap[String,Int](numAtoms);

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

            val globalIndex = Long.parseLong(line.substring(15n,20n).trim());
            // multiply coords by 10 to convert nm to Angstroms
            val center = Point3d(Double.parseDouble(line.substring(20n,28n)),
                                 Double.parseDouble(line.substring(28n,36n)),
                                 Double.parseDouble(line.substring(36n,44n)));
            particleData.addAtom(globalIndex, atomTypeIndex, center);

            if (line.length() >= 67) {
                // velocities are included
                particleData.dx(i) = Vector3d(Double.parseDouble(line.substring(44n,52n)),
                                              Double.parseDouble(line.substring(52n,60n)),
                                              Double.parseDouble(line.substring(60n,68n)));
            }
            val particleKey = residueName+residueNumber+atomName;
            particleMap.put(particleKey, i);
        }

        // read box edge lengths
        val boxLine = file.readLine();
        val x = Double.parseDouble(boxLine.substring( 0n,10n).trim());
        val y = Double.parseDouble(boxLine.substring(10n,20n).trim());
        val z = Double.parseDouble(boxLine.substring(20n,30n).trim());
        particleData.boxEdges = Vector3d(x,y,z);
        
        file.close();

        createBonds(particleData, particleMap);
    }

    /** 
     * Create bonds for each residue from the bond structure held in
     * particleData.moleculeTypes
     */
    private def createBonds(particleData:AnummParticleData, atomMap:HashMap[String,Int]) {
        val bonds = new ArrayList[Bond]();
        val angles = new ArrayList[BondAngle]();
        for (i in 0..(particleData.numAtoms()-1)) {
            val atomTypeIndex = particleData.atomTypeIndex(i);
            val residueNumber = particleData.residueNumber(i);
            val moleculeType = particleData.moleculeTypes(particleData.moleculeTypeIndex(residueNumber));
            if (moleculeType.bonds != null) {
                for (bond in moleculeType.bonds) {
                    // bonds follow first atom
                    if (atomTypeIndex == moleculeType.atomTypes(bond.atom1Index-1)) {
                        val atom1Key = moleculeType.name+residueNumber+moleculeType.atomNames(bond.atom1Index-1);
                        val atom1Index = atomMap.get(atom1Key);
                        val atom2Key = moleculeType.name+residueNumber+moleculeType.atomNames(bond.atom2Index-1);
                        val atom2Index = atomMap.get(atom2Key);
                        //Console.OUT.println("creating bond from " + atom1Key + " index " + atom1Index + " to " + atom2Key + " index " + atom2Index + " type " + bond.typeIndex);
                        bonds.add(Bond(atom1Index, atom2Index, bond.typeIndex));
                    }
                }
            }
            if (moleculeType.angles != null) {
                for (angle in moleculeType.angles) {
                    // bond angles follow second atom
                    if (atomTypeIndex == moleculeType.atomTypes(angle.atom2Index-1)) {
                        val atom1Key = moleculeType.name+residueNumber+moleculeType.atomNames(angle.atom1Index-1);
                        val atom1Index = atomMap.get(atom1Key);
                        val atom2Key = moleculeType.name+residueNumber+moleculeType.atomNames(angle.atom2Index-1);
                        val atom2Index = atomMap.get(atom2Key);
                        val atom3Key = moleculeType.name+residueNumber+moleculeType.atomNames(angle.atom3Index-1);
                        val atom3Index = atomMap.get(atom3Key);
                        //Console.OUT.println("creating angle " + atom1Key + " - " + atom2Key + " - " + atom3Key + " type " + angle.typeIndex);
                        angles.add(BondAngle(atom1Index, atom2Index, atom3Index, angle.typeIndex));
                    }
                }
            }
        }
        particleData.bonds = bonds.toRail();
        particleData.angles = angles.toRail();
    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

