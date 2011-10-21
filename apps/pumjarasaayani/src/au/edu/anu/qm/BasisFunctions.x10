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
import au.edu.anu.chem.Molecule;

/**
 * BasisFunctions.x10
 *
 * Generates basis functions for a given Molecule object
 *
 * @author: V.Ganesh
 */
public struct BasisFunctions { 
    public val molecule:Molecule[QMAtom];
    public val basisName:String;
    public val basisFunctions:ArrayList[ContractedGaussian];
    public val shellList:ShellList;

    public def this(mol:Molecule[QMAtom], basNam:String, basisDir:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian](); 
        initBasisFunctions(basisDir);

        shellList = new ShellList(mol);
    }

    private def initBasisFunctions(basisDir:String) {
        val basisSet:BasisSet = new BasisSet(basisName, basisDir);
        var intIndx:Int = 0;
        val plInst = PowerList.getInstance();

        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom      = molecule.getAtom(atmno);
            val atomBasis = basisSet.getBasis(atom);
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
                    val cg = new ContractedGaussian(centre, power, primitives, intIndx, true);
                    basisFunctions.add(cg);
                } // end for

                val am = orb.angularMomentum;
                val atomPrimitives = new Array[PrimitiveGaussian](coeffs.size);
                for(var i:Int=0; i<coeffs.size; i++) {
                    atomPrimitives(i) = new PrimitiveGaussian(centre, Power(am, 0, 0), exps(i), coeffs(i), false);
                }
                val acg = new ContractedGaussian(centre, Power(am, 0, 0), atomPrimitives, intIndx, false);
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
                var lm:Int = bfi.getMaximumAngularMomentum();
                var denom:Double = 1.0;

                while(lm>1) { denom *= (2.0*lm-1.0); lm--; }
              
                val lmn   = bfi.getMaximumAngularMomentum();
                val fact2 = Math.pow(2.0, lmn) / Math.pow(denom, 0.5);
              
                val prms  = bfi.getPrimitives(); 
                for(var j:Int=0; j<prms.size; j++) {
                    val prmj = prms(j); 
                    val fact3 = Math.pow(prmj.exponent, (2.0*lmn+3.0)/4.0);
                    prmj.setCoefficient(prmj.coefficient*fact1*fact2*fact3);
                } // end for

                val pi3O2By2Tlmn = Math.pow(Math.PI, 1.5) / Math.pow(2.0, lmn); 
                val factC = denom * pi3O2By2Tlmn;
                var factA:Double = 0.0, factB:Double = 0.0;
                var coefExpoSum:Double = 0.0;

                for(var k:Int=0; k<prms.size; k++) {
                   val prmk = prms(k);
                   for(var l:Int=0; l<prms.size; l++) {
                       val prml = prms(l);

                       factA = prmk.coefficient * prml.coefficient;
                       factB = Math.pow(prmk.exponent + prml.exponent, lmn+1.5);
                       coefExpoSum += factA/factB;
                   } // end for
                } // end for 

                val norm = Math.pow(coefExpoSum*factC, -0.5);
                for(var j:Int=0; j<prms.size; j++) {
                    val prmj = prms(j); 
                    prmj.setCoefficient(prmj.coefficient*norm);
                } // end for 
            } // end for
        } // end for
    }

    public def getBasisName() = this.basisName;
    public def getBasisFunctions() = basisFunctions;
    public def getShellList() = shellList;
}
