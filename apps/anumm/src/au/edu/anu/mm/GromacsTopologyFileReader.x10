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

import au.edu.anu.chem.mm.AtomType;

/**
 * This class reads GROMACS topology files (*.top) to populate particle
 * data and force field parameters.
 * @see http://manual.gromacs.org/current/online/gro.html
 */
public class GromacsTopologyFileReader { 
    var fileName:String;
    var currentLine:String;

    val atomTypes = new ArrayList[AtomType]();
    val moleculeTypes = new ArrayList[MoleculeType]();
    val bondTypes = new ArrayList[BondStretchParameters]();
    val angleTypes = new ArrayList[BondAngleParameters]();

    public def this(fileName:String) { 
        this.fileName = fileName;
    }

    public def readMoleculeType(file:Reader, allAtomTypes:Rail[AtomType]) {
        currentLine = file.readLine();
        while (!finishedSection()) {
            if (!isComment()) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val moleculeType = new MoleculeType();
                moleculeType.name = words(0);

                while (currentLine != null) {
                    if (currentLine.indexOf("[ atoms ]") >= 0) {
                        readAtoms(file, moleculeType, allAtomTypes);
                    } else if (currentLine.indexOf("[ bonds ]") >= 0) {
                        readBonds(file, moleculeType);
                    } else if (currentLine.indexOf("[ pairs ]") >= 0) {
                        skipSection(file);
                    } else if (currentLine.indexOf("[ angles ]") >= 0) {
                        readAngles(file, moleculeType);
                    } else if (currentLine.indexOf("[ dihedrals ]") >= 0) {
                        skipSection(file);
                    } else if (currentLine.indexOf("[") >= 0) {
                        // finished this molecule type
                        break;
                    } else {
                        currentLine = file.readLine();
                    }
                }
                moleculeTypes.add(moleculeType);
            } else {
                currentLine = file.readLine();
            }
        }
    }

    public def readAtoms(file:Reader, molecule:MoleculeType, allAtomTypes:Rail[AtomType]) {
        currentLine = file.readLine();
        molecule.atomNames = new ArrayList[String]();
        molecule.atomTypes = new ArrayList[Int]();
        while (!finishedSection()) {
            if (!isComment() && currentLine.trim().length() > 0n) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val atomNumber = Int.parseInt(words(0));
                val atomTypeName = words(1);
                val atomName = words(4);
                val charge = Double.parseDouble(words(6));
                var mass:Double = 0.0;
                var atomicNumber:Int = -1n;
                for (j in 0..(allAtomTypes.size-1)) {
                    val atomType = allAtomTypes(j);
                    if (atomTypeName.equals(atomType.name)) {
                        mass = atomType.mass;
                        atomicNumber = atomType.atomicNumber;
                        //Console.OUT.println("found atomType " + atomTypeName + " " + mass);
                        atomTypes.add(AtomType(atomTypeName, atomicNumber, mass, charge));

                        // TODO link atomName and index in MoleculeType
                        molecule.atomNames.add(atomName);
                        molecule.atomTypes.add((atomTypes.size()-1) as Int);
                        break;
                    }
                }
                if (atomicNumber == -1n) {
                    Console.ERR.println("did not find atom type " + atomTypeName);
                }
            }
            currentLine = file.readLine();
        }
    }

    public def readBonds(file:Reader, molecule:MoleculeType) {
        val bonds = new ArrayList[Bond]();
        currentLine = file.readLine();
        while (!finishedSection()) {
            if (!isComment() && currentLine.trim().length() > 0n) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val atom1Index = Int.parseInt(words(0));
                val atom2Index = Int.parseInt(words(1));
                val functionType = Int.parseInt(words(2)); // not used

                val idealRadius = Double.parseDouble(words(3));
                val forceConstant = Double.parseDouble(words(4));
                val atom1TypeIndex = molecule.atomTypes(atom1Index-1);
                val atom2TypeIndex = molecule.atomTypes(atom2Index-1);
                bondTypes.add(new BondStretchParameters(atom1TypeIndex, atom2TypeIndex, idealRadius, forceConstant));
                //Console.OUT.println("adding new bond stretch params " + atom1TypeIndex + " to " + atom2TypeIndex + " idealRadius " + idealRadius);

                //Console.OUT.println("molecule " + molecule.name + " has bond " + atom1Index + " to " + atom2Index + " type " + (bondTypes.size()-1));
                bonds.add(new Bond(atom1Index, atom2Index, (bondTypes.size()-1) as Int));
            }

            currentLine = file.readLine();
        }
        molecule.bonds = bonds;
    }

    public def readAngles(file:Reader, molecule:MoleculeType) {
        val angles = new ArrayList[BondAngle]();
        currentLine = file.readLine();
        while (!finishedSection()) {
            if (!isComment() && currentLine.trim().length() > 0n) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val atom1Index = Int.parseInt(words(0));
                val atom2Index = Int.parseInt(words(1));
                val atom3Index = Int.parseInt(words(2));
                val function = Int.parseInt(words(3)); // not used
                 // convert to radians
                val idealAngle = Double.parseDouble(words(4)) / 180.0 * Math.PI;
                val forceConstant = Double.parseDouble(words(5));

                //Console.OUT.println("molecule " + molecule.name + " has angle " + atom1Index + " - " + atom2Index + " - " + atom3Index + " idealAngle " + idealAngle + " forceConstant " + forceConstant);
                val atom1TypeIndex = molecule.atomTypes(atom1Index-1);
                val atom2TypeIndex = molecule.atomTypes(atom2Index-1);
                val atom3TypeIndex = molecule.atomTypes(atom3Index-1);
                angleTypes.add(new BondAngleParameters(atom1TypeIndex, atom2TypeIndex, atom3TypeIndex, idealAngle, forceConstant));

                angles.add(new BondAngle(atom1Index, atom2Index, atom3Index, (angleTypes.size()-1) as Int));
            }

            currentLine = file.readLine();
        }
        molecule.angles = angles;
    }

    public def readAtomTypes(file:Reader):Rail[AtomType] {
        val atomTypes = new ArrayList[AtomType]();
        currentLine = file.readLine();
        while (!finishedSection()) {
            if (!isComment() && currentLine.trim().length() > 0n) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val name = words(0);
                val atomicNumber = Int.parseInt(words(1));
                val mass = Double.parseDouble(words(2));
                val charge = Double.parseDouble(words(3));
                atomTypes.add(new AtomType(name, atomicNumber, mass, charge));
            }

            currentLine = file.readLine();
        }
        return atomTypes.toRail();
    }

    public def readMolecules(file:Reader, moleculeTypeIndex:ArrayList[Int]) {
        currentLine = file.readLine();
        while (!finishedSection()) {
            if (!isComment()) {
                val words = StringSplitter.splitOnWhitespace(currentLine);
                val moleculeTypeName = words(0);
                val number = Long.parseLong(words(1));
                moleculeTypeIndex.resize((moleculeTypeIndex.size() + number), -1n);
            }
            currentLine = file.readLine();
        }
    }

    public def readTopology(particleData:AnummParticleData) {
        // perform pre-processing
        val gmxLib = System.getenv("GMXLIB");
        if (gmxLib == null) {
            throw new UnsupportedOperationException("GromacsTopologyFileReader: GMXLIB is not set");
        }
        // construct GROMACS topology file using C preprocessor,
        // inhibiting generation of linemarkers
        val processedFile = Runtime.execForRead("cpp -P -I"+gmxLib+" "+fileName);

        var allAtomTypes:Rail[AtomType] = null;
        // add a dummy type as element 0, because GROMACS uses 1-based indexing
        val dummy = new MoleculeType();
        dummy.name = "unknown";
        moleculeTypes.add(dummy);
        val moleculeTypeIndex = new ArrayList[Int]();
        moleculeTypeIndex.add(-1n);

        currentLine = processedFile.readLine();
        try {
            while (currentLine != null) {
                if (!isComment()) {
                    if (currentLine.indexOf("[ atomtypes ]") >= 0) {
                        allAtomTypes = readAtomTypes(processedFile);
                    } else if (currentLine.indexOf("[ moleculetype ]") >= 0) {
                        readMoleculeType(processedFile, allAtomTypes);
                    } else if (currentLine.indexOf("[ molecules ]") >= 0) {
                        readMolecules(processedFile, moleculeTypeIndex);
                    } else {
                        currentLine = processedFile.readLine();
                    }
                } else {
                    currentLine = processedFile.readLine();
                }
            }
        } catch (e:EOFException) {
            // no more lines
        }

/*
        Console.OUT.println("found simulation atom types:");
        for (i in 0..(atomTypes.size()-1)) {
            Console.OUT.println(atomTypes(i));
        }
*/
        particleData.atomTypes = atomTypes.toRail();
        particleData.moleculeTypes = moleculeTypes.toRail();
        particleData.moleculeTypeIndex = moleculeTypeIndex.toRail();

        particleData.bondTypes = bondTypes.toRail();
        particleData.angleTypes = angleTypes.toRail();

        processedFile.close();
    }

    private def isComment() {
        return currentLine.startsWith(";");
    }

    private def finishedSection() {
        if (currentLine == null) return true;
        else {
            val trimmed = currentLine.trim();
            if (trimmed.startsWith("["))
                return true;
        }
        return false;
    }

    private def skipSection(file:Reader) {
        currentLine = file.readLine();
        while (!finishedSection()) {
            currentLine = file.readLine();
        }
    }

    public def setFileName(fileName:String) {
        this.fileName = fileName;
    }
}

