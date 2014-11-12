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

import au.edu.anu.util.StringSplitter;

import x10.compiler.Inline;
import x10.io.File;
import x10.io.Reader;
import x10.io.EOFException;
import x10.util.ArrayList;

/**
 * This class reads GROMACS topology files (*.top) of the following format:
 *  <title>
 *  <number of atoms>
 *  <residue number> <residue name> <atom name> <atom number> x y z [vx vy vz]// repeated
 *  <box x> <box y> <box z>
 * TODO read force field data from processed topology
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsTopologyFileReader { 
    var fileName:String;

    public def this(fileName:String) { 
        this.fileName = fileName;
    }

    public def readResidue(file:Reader, particleData:ParticleData, forceField:ForceField, residueType:ArrayList[String]) {
        var line:String = file.readLine();
        while (line != null && line.trim().length() > 0n) {
            if (!isComment(line)) {
                val words = StringSplitter.splitOnWhitespace(line);
                val residueName = words(0);
                residueType.add(residueName);

                break; // TODO read bond stretch, angle parameters
            }
            line = file.readLine();
        }
    }

    public def readMolecules(file:Reader, particleData:ParticleData, forceField:ForceField, residueTypeIndex:ArrayList[Int]) {
        var line:String = file.readLine();
        if (isComment(line) || line.trim().length() == 0n) line = file.readLine(); // may be a comment
        while (line != null && line.trim().length() > 0n) {
            if (!isComment(line)) {
                val words = StringSplitter.splitOnWhitespace(line);
                val residueName = words(0);
                val number = Long.parseLong(words(1));
                residueTypeIndex.resize((residueTypeIndex.size() + number), -1n);
            }
            line = file.readLine();
        }
    }

    public def readBonds(file:Reader, particleData:ParticleData, forceField:ForceField) {
        val bonds = new ArrayList[Bond]();
/*   
        var line:String = file.readLine();
        var i:Long = 0;
        while (line != null && line.trim().length() > 0n) {
            if (!isComment(line)) {

                val words = StringSplitter.splitOnWhitespace(line);
                val atom1Index = Int.parseInt(words(0));//Int.parseInt(line.substring(0n,3n).trim())-1;
                val atom2Index = Int.parseInt(words(1));//Int.parseInt(line.substring(3n,6n).trim())-1;
                val bondType = Int.parseInt(words(2));//Int.parseInt(line.substring(6n,9n).trim());
                Console.OUT.println("bond " + atom1Index + " to " + atom2Index + " type " + bondType);
                bonds.add(new Bond(atom1Index, atom2Index, forceField.getBondTypeIndex(particleData.species(atom1Index), particleData.species(atom2Index), bondType)));

            }

            line = file.readLine();
        }
*/
        particleData.bonds = bonds.toRail();
    }

    public def readTopology(particleData:ParticleData, forceField:ForceField) {
        // perform pre-processing
        val gmxLib = System.getenv("GMXLIB");
        if (gmxLib == null) {
            throw new UnsupportedOperationException("GromacsTopologyFileReader: GMXLIB is not set");
        }
        val processedFile = Runtime.execForRead("cpp -I"+gmxLib+" "+fileName);

        val residueType = new ArrayList[String](); 
        val residueTypeIndex = new ArrayList[Int](); 

        var line:String = processedFile.readLine();
        try {
            while (line != null) {
                if (!isComment(line)) {
                    if (line.indexOf("[ bonds ]") >= 0) {
                        readBonds(processedFile, particleData, forceField);
                    } else if (line.indexOf("[ molecules ]") >= 0) {
                        readMolecules(processedFile, particleData, forceField, residueTypeIndex);
                    } else if (line.indexOf("[ moleculetype ]") >= 0) {
                        readResidue(processedFile, particleData, forceField, residueType);
                    }
                }
                line = processedFile.readLine();
            }
        } catch (e:EOFException) {
            // no more lines
        }
/*
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
*/
        particleData.residueType = residueType.toRail();
        particleData.residueTypeIndex = residueTypeIndex.toRail();
        processedFile.close();
    }

    private @Inline def isComment(line:String) {
        return line.startsWith(";");
    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

