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

    public def this(name:String) { 
       this.name = name;

       basisInfo = new HashMap[String, AtomicBasis]();
       init(name);
    } 

    public def this(name:String, basisDir:String) {
       Console.OUT.println("\tReading in basis info. for " + name + " from " + basisDir);

       basisInfo = new HashMap[String, AtomicBasis]();
       this.name = name;

       try {
         init(name, basisDir);
       } catch(e:Exception) {
         Console.OUT.println("Unable to read basis from : " + basisDir + ".");
         Console.OUT.println("Will use sto3g basis.");

         init("sto3g");
       } // end of try .. catch block
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
 
    private def init(name:String) {
        if (name.equals("sto3g")) {

            var orbType:String = "S";
            var exps:Rail[Double] = [3.425251, 0.623914, 0.168855];
            var coeffs:Rail[Double] = [0.154329, 0.535328, 0.444635];
            val hOrb = new Orbital(orbType, exps, coeffs);
            var orbs : Rail[Orbital] = [hOrb as Orbital];
            val hBasis = new AtomicBasis(orbs);

            basisInfo.put("H", hBasis);

            orbType = "S";
            exps = [71.616837, 13.045096, 3.530512];
            coeffs = [0.154329, 0.535328, 0.444635];
            val cOrb1 = new Orbital(orbType, exps, coeffs);
            orbType = "S";
            exps = [2.941249, 0.683483, 0.222290];
            coeffs = [-0.099967, 0.399513, 0.700115];
            val cOrb2 = new Orbital(orbType, exps, coeffs);
            orbType = "P";
            exps = [2.941249, 0.683483, 0.222290];
            coeffs = [0.155916, 0.607684, 0.391957];
            val cOrb3 = new Orbital(orbType, exps, coeffs);
            orbs = [cOrb1, cOrb2, cOrb3];
            val cBasis = new AtomicBasis(orbs);

            basisInfo.put("C", cBasis);

            orbType = "S";
            exps = [99.106169, 18.052312, 4.885660];
            coeffs = [0.154329, 0.535328, 0.444635];
            val nOrb1 = new Orbital(orbType, exps, coeffs);
            orbType = "S";
            exps = [3.780456, 0.878497, 0.285714];
            coeffs = [-0.099967, 0.399513, 0.700115];
            val nOrb2 = new Orbital(orbType, exps, coeffs);
            orbType = "P";
            exps = [3.780456, 0.878497, 0.285714];
            coeffs = [0.155916, 0.607684, 0.391957];
            val nOrb3 = new Orbital(orbType, exps, coeffs);
            orbs = [nOrb1, nOrb2, nOrb3];
            val nBasis = new AtomicBasis(orbs);

            basisInfo.put("N", nBasis);

            orbType = "S";
            exps = [130.709321, 23.808866, 6.443608];
            coeffs = [0.154329, 0.535328, 0.444635];
            val oOrb1 = new Orbital(orbType, exps, coeffs);
            orbType = "S";
            exps = [5.033151, 1.169596, 0.380389];
            coeffs = [-0.099967, 0.399513, 0.700115];
            val oOrb2 = new Orbital(orbType, exps, coeffs);
            orbType = "P";
            exps = [5.033151, 1.169596, 0.380389];
            coeffs = [0.155916, 0.607684, 0.391957];
            val oOrb3 = new Orbital(orbType, exps, coeffs);
            orbs = [oOrb1, oOrb2, oOrb3];
            val oBasis = new AtomicBasis(orbs);

            basisInfo.put("O", oBasis);    
       } // end if
    }

    public def getName() = this.name;

    public def getBasis(atom:Atom) = basisInfo.getOrElse(atom.symbol, null);
}

