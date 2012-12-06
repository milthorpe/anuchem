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

import x10.util.ArrayList;

import x10.matrix.DenseMatrix;
import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

/**
 * Generates basis functions for a given Molecule object
 *
 * @author: V.Ganesh
 */
public struct BasisFunctions { 
    public val molecule:Molecule[QMAtom];
    public val basisName:String;
    public val basisFunctions:ArrayList[ContractedGaussian];
    public val shellList:ShellList;
    public val SADMatrix:DenseMatrix;

    public def this(mol:Molecule[QMAtom], basNam:String, basisDir:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian](); 
        val basisSet:BasisSet = new BasisSet(basisName, basisDir);

        val size = initBasisFunctions(basisSet);
        shellList = new ShellList(mol);
        SADMatrix = new DenseMatrix(size,size);
        val jd = JobDefaults.getInstance();
        if (jd.guess.equals(JobDefaults.GUESS_SAD))     
            initDensity(basisSet);
    }

    private def initDensity(basisSet:BasisSet) {
        var shift:Int=0; 
        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom = molecule.getAtom(atmno);
            val aDensity = basisSet.getDensity(atom);
            if (aDensity == null) {
                throw new Exception("No density matrix found for atom type " + atom.species);
            }
            val matsize = aDensity.M;
            for ([i,j] in 0..(matsize-1)*0..(matsize-1)) {
                SADMatrix(i+shift,j+shift) = aDensity(i,j);
            }
            shift+=matsize;
        }
    }

    private def initBasisFunctions(basisSet:BasisSet):Int {
        
        var intIndx:Int = 0;
        val plInst = PowerList.getInstance();

        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom      = molecule.getAtom(atmno);
            val atomBasis = basisSet.getBasis(atom);
            if (atomBasis == null) {
                throw new Exception("No basis found for atom type " + atom.species);
            }
            val orbitals  = atomBasis.orbitals;
            val atombfs   = new ArrayList[ContractedGaussian]();

            for(var orbno:Int=0; orbno<orbitals.size; orbno++) { 
                val orb = orbitals(orbno);
                val shape = orb.shape;
                val pList = plInst.getPowers(shape);

                val coeffs = orb.coefficients;
                val exps = orb.exponents;
                val centre = atom.centre;

                for(var l:Int=0; l<pList.size; l++) {
                    val power = pList(l);
                    val primitives = new Array[PrimitiveGaussian](coeffs.size);
                    for(var i:Int=0; i<coeffs.size; i++) {
                        primitives(i) = new PrimitiveGaussian(centre, power, exps(i), coeffs(i), true);
                    }
                    val cg = new ContractedGaussian(centre, power, exps, coeffs, intIndx, true);
                    basisFunctions.add(cg);
                } // end for

                val am = orb.angularMomentum;
                val atomPrimitives = new Array[PrimitiveGaussian](coeffs.size);
                for(var i:Int=0; i<coeffs.size; i++) {
                    atomPrimitives(i) = new PrimitiveGaussian(centre, Power(am, 0, 0), exps(i), coeffs(i), false);
                }
                val atomCoeffs = new Array[Double](coeffs); // note: atom coefficents are subsequently normalized
                val acg = new ContractedGaussian(centre, Power(am, 0, 0), exps, atomCoeffs, intIndx, false);
                atombfs.add(acg);

                intIndx += ((am+1)*(am+2)/2);
            } // end for

            atom.setBasisFunctions(atombfs);  
        } // end for

        // normalization of atom basis functions : mostly from Alistair's code, this essentially same as the above 
        // code, except is performed over a differently arranged set of ContractedGaussian's 
        val fact1 = Math.pow(2.0/Math.PI, 0.75);
        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom = molecule.getAtom(atmno);
            val bfs  = atom.getBasisFunctions();
            val nbf  = bfs.size();

            for(var i:Int=0; i<nbf; i++) {
                val bfi    = bfs.get(i);
                val lmn   = bfi.getMaximumAngularMomentum();
                var lm:Int = lmn;
                var denom:Double = 1.0;

                while(lm>1) { denom *= (2.0*lm-1.0); lm--; }
              
                val fact2 = Math.pow2(lmn) / Math.pow(denom, 0.5);
              
                val exponents = bfi.exponents;
                val coefficients = bfi.coefficients;
                for(var j:Int=0; j<coefficients.size; j++) {
                    val fact3 = Math.pow(exponents(j), (2.0*lmn+3.0)/4.0);
                    coefficients(j) = coefficients(j)*fact1*fact2*fact3;
                } // end for

                val pi3O2By2Tlmn = Math.pow(Math.PI, 1.5) / Math.pow2(lmn);
                val factC = denom * pi3O2By2Tlmn;
                var factA:Double = 0.0, factB:Double = 0.0;
                var coefExpoSum:Double = 0.0;

                for(var k:Int=0; k<coefficients.size; k++) {
                   for(var l:Int=0; l<coefficients.size; l++) {
                       factA = coefficients(k) * coefficients(l);
                       factB = Math.pow(exponents(k) + exponents(l), lmn+1.5);
                       coefExpoSum += factA/factB;
                   } // end for
                } // end for 

                val norm = Math.pow(coefExpoSum*factC, -0.5);
                for(var j:Int=0; j<coefficients.size; j++) {
                    coefficients(j) *= norm;
                } // end for 
            } // end for
        } // end for	
	    return intIndx;
    }

    public def getNormalizationFactors():Rail[Double] {
        val norms = new Array[Double](basisFunctions.size());
        for(var i:Int=0; i<basisFunctions.size(); i++) {
            val power = basisFunctions(i).power;
             for(var j:Int=0; j<1; j++) {
                var lmx:Int = power.getL();
                var lmy:Int = power.getM();
                var lmz:Int = power.getN();
                var lmt:Int = power.getTotalAngularMomentum();
                var denom:Double=1.0;

                while(lmt>1) { denom *= (2.0*lmt-1.0); lmt--; }
                while(lmx>1) { denom /= (2.0*lmx-1.0); lmx--; }
                while(lmy>1) { denom /= (2.0*lmy-1.0); lmy--; }
                while(lmz>1) { denom /= (2.0*lmz-1.0); lmz--; }

                norms(i) = Math.sqrt(denom);
            }
        }
        return norms;
    }

    public def getBasisName() = this.basisName;
    public def getBasisFunctions() = basisFunctions;
    /** @return the guess for the density matrix formed from a superposition of atomic densities */
    public def getSAD() = SADMatrix;
    public def getShellList() = shellList;
}

