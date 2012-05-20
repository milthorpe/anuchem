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

import x10.io.EOFException;
import x10.io.File;
import x10.io.FileReader;
import x10.util.ArrayList;
import x10.util.HashMap;
import au.edu.anu.chem.Atom;
import au.edu.anu.util.StringSplitter;
import x10.matrix.DenseMatrix;

/**
 * Represents a basis set, used for setting up basis functions (by BasisFunctions)
 * Basis set input files are expected to be in Gaussian 94 format.
 * @see http://www.emsl.pnl.gov/forms/basisform.html
 *
 * @author V. Ganesh, milthorpe
 */
public class BasisSet { 
    val name:String;
    val basisInfo:HashMap[String, AtomicBasis];
    val basisAtomicDensity:HashMap[String, DenseMatrix];
    val roZ:Double;

    public def this(name:String, basisDir:String) {
        Console.OUT.println("\tReading in basis info. for " + name + " from " + basisDir);

        basisInfo = new HashMap[String, AtomicBasis]();
	    basisAtomicDensity = new HashMap[String, DenseMatrix]();
 
        this.name = name;
        val jd = JobDefaults.getInstance();
        this.roZ=jd.roZ;
        try {
            init(name, basisDir);
        } catch(e:Exception) {
            throw new Exception("Unable to read basis from : "+basisDir, e);
        }

        if (jd.guess.equals(JobDefaults.GUESS_SAD)) {
            try {
                initDensity(name, basisDir);
            } catch(e:Exception) {
                throw new Exception("Unable to read density from : "+basisDir, e);
            } 
        }
    }

    /** 
     * Initialize this basis set from a file in Gaussian 94 format. 
     * Limitations:
     * - only supports atomic symbol specifications (not center numbers)
     * - assumes only one shell block per atom type
     * @see http://www.gaussian.com/g_tech/g_ur/k_gen.htm
     */
    private def init(name:String, basisDir:String) {
        val fileName = basisDir + File.SEPARATOR + name;
        val fil = new FileReader(new File(fileName));       

        while(readCenterDefinitionBlock(fil) != null);

        fil.close();
    }

    private def isComment(line:String) {
        return line.startsWith("!");
    }

    private def isSectionDivider(line:String) {
        return line.startsWith("****") || line.startsWith("++++");
    }

    /** Read the center definition block for a single atom type. */
    private def readCenterDefinitionBlock(fil:FileReader):AtomicBasis {
        try {
            var line:String = fil.readLine();
            while (line != null) {
                line = line.trim();
                if (line.length() > 0 && !(isComment(line) || isSectionDivider(line))) {
                    val centerIdentiferWords = StringSplitter.splitOnWhitespace(line);
                    val symbol = centerIdentiferWords(0);

                    val orbitalList = new ArrayList[Orbital]();
                    line = fil.readLine();
                    while(!isSectionDivider(line)) {
                        val shellWords = StringSplitter.splitOnWhitespace(line);
                        val shellType = shellWords(0);
                        val numGaussians = Int.parseInt(shellWords(1));
                        val scaleFactor = Double.parseDouble(shellWords(2)); // TODO what to do with scaleFactor

                        if (shellType.equals("SP")) {
                            val exps = new Array[Double](numGaussians);
                            val sCoeffs = new Array[Double](numGaussians);
                            val pCoeffs = new Array[Double](numGaussians);
                            for (i in 0..(numGaussians-1)) {
                                line = fil.readLine();
                                val gaussianWords = StringSplitter.splitOnWhitespace(line);
                                exps(i) = Double.parseDouble(gaussianWords(0))*roZ*roZ;
                                sCoeffs(i) = Double.parseDouble(gaussianWords(1));
                                pCoeffs(i) = Double.parseDouble(gaussianWords(2));
                            }
                            orbitalList.add(new Orbital("S", exps, sCoeffs));
                            orbitalList.add(new Orbital("P", exps, pCoeffs));
                        } else {
                            val exps = new Array[Double](numGaussians);
                            val coeffs = new Array[Double](numGaussians);
                            for (i in 0..(numGaussians-1)) {
                                line = fil.readLine();
                                val gaussianWords = StringSplitter.splitOnWhitespace(line);
                                exps(i) = Double.parseDouble(gaussianWords(0))*roZ*roZ;
                                coeffs(i) = Double.parseDouble(gaussianWords(1));
                            }
                            orbitalList.add(new Orbital(shellType, exps, coeffs));
                        }
                        line = fil.readLine();
                    }

                    val atomBasis = new AtomicBasis(orbitalList.toArray());

                    basisInfo.put(symbol, atomBasis);

                    return atomBasis;
                }
                line = fil.readLine();
            }
        } catch (e:EOFException) {
            // no more center definition blocks
        }

        return null;
    }

    private def initDensity(name:String, basisDir:String) {
        val fileName = basisDir + File.SEPARATOR + name + ".P";
        val fil = new FileReader(new File(fileName));       
        val noOfAtoms = Int.parseInt(fil.readLine());

        for(var i:Int=0; i<noOfAtoms; i++) {
            val words = StringSplitter.splitOnWhitespace(fil.readLine());
            val symbol = words(0);
            val noOfFunctions = Int.parseInt(words(1));

            val density = new DenseMatrix(noOfFunctions, noOfFunctions);
            for(var j:Int=0; j<noOfFunctions; j++) {
                val words1 = StringSplitter.splitOnWhitespace(fil.readLine());
                for(var k:Int=0; k<noOfFunctions; k++) {
                    density(j,k) = Double.parseDouble(words1(k));
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

