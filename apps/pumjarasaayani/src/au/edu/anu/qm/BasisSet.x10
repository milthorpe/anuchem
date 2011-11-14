/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.io.File;
import x10.io.FileReader;
import x10.util.HashMap;
import au.edu.anu.chem.Atom;
import x10x.matrix.Matrix;

/**
 * BasisSet.x10
 *
 * Represents a basis set, used for setting up basis functions (by BasisFunctions)
 *
 * @author: V.Ganesh
 */
public class BasisSet { 
    val name:String;
    val basisInfo:HashMap[String, AtomicBasis];
    val basisAtomicDensity:HashMap[String, Matrix];

    public def this(name:String, basisDir:String) {
        Console.OUT.println("\tReading in basis info. for " + name + " from " + basisDir);

        basisInfo = new HashMap[String, AtomicBasis]();
	    basisAtomicDensity = new HashMap[String, Matrix]();
 
        this.name = name;
        try {
            init(name, basisDir);
        } catch(e:Exception) {
            throw new Exception("Unable to read basis from : "+basisDir, e);
        }

        try {
            initDensity(name, basisDir);
        } catch(e:Exception) {
            throw new Exception("Unable to read density from : "+basisDir, e);
        }
    }

    private def init(name:String, basisDir:String) {
        val fileName = basisDir + File.SEPARATOR + name;
        val fil = new FileReader(new File(fileName));       
        val noOfAtoms = Int.parseInt(fil.readLine());

        for(var i:Int=0; i<noOfAtoms; i++) {
            val words = fil.readLine().split(" ");
            val symbol = words(0);
            val noOfContractions = Int.parseInt(words(1));

            val orbitals = new Array[Orbital](noOfContractions);
            for(var j:Int=0; j<noOfContractions; j++) {
                val words1 = fil.readLine().split(" ");
                var orbitalType:String = words1(0);
                var noOfPrimitives:Int = Int.parseInt(words1(1));

                val exps:Rail[Double] = new Array[Double](noOfPrimitives);
                val coeffs:Rail[Double] = new Array[Double](noOfPrimitives);
                for(var k:Int=0; k<noOfPrimitives; k++) { 
                    val words2 = fil.readLine().split(" ");
                    exps(k) = Double.parseDouble(words2(0));
                    coeffs(k) = Double.parseDouble(words2(1));
                }
                orbitals(j) = new Orbital(orbitalType, exps, coeffs); 
            }

            val atomBasis = new AtomicBasis(orbitals);

            basisInfo.put(symbol, atomBasis);
        } // end for

        fil.close();
    }

    private def initDensity(name:String, basisDir:String) {
        val fileName = basisDir + File.SEPARATOR + name + ".P";
        val fil = new FileReader(new File(fileName));       
        val noOfAtoms = Int.parseInt(fil.readLine());

        for(var i:Int=0; i<noOfAtoms; i++) {
            val words = fil.readLine().split(" ");
            val symbol = words(0);
            val noOfFunctions = Int.parseInt(words(1));

            val density = new Matrix(noOfFunctions, noOfFunctions);
	        val aDensity = density.getMatrix();
            for(var j:Int=0; j<noOfFunctions; j++) {
                val words1 = fil.readLine().split(" ");
                for(var k:Int=0; k<noOfFunctions; k++) {
                    aDensity(j,k) = Double.parseDouble(words1(k));
                }
            }

            basisAtomicDensity.put(symbol, density);
        } // end for

        fil.close();
    }

    public def getName() = this.name;

    public def getBasis(atom:Atom) = basisInfo.getOrElse(atom.symbol, null);
    public def getDensity(atom:Atom) = basisAtomicDensity.getOrElse(atom.symbol, null);

}

