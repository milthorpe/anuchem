/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

import x10.util.ArrayList;
import au.edu.anu.chem.Molecule;

/**
 * BasisFunctions.x10
 *
 * Generates basis functions for a given Molecule object
 *
 * @author: V.Ganesh
 */
public class BasisFunctions { 
    global val molecule:Molecule[QMAtom]!;
    global val basisName:String;
    global val basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    global val shellList:ShellList!;

    public def this(mol:Molecule[QMAtom]!, basNam:String, basisDir:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian{self.at(this)}](); 
        initBasisFunctions(basisDir);

        shellList = new ShellList();
        initShellList();
    } 

    private def initBasisFunctions(basisDir:String) {
        val basisSet:BasisSet{self.at(this)} = new BasisSet(basisName, basisDir);
        var indx:Int = 0;
        var intIndx:Int = 0;
        val plInst = PowerList.getInstance();

        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom      = molecule.getAtom(atmno);
            val atomBasis = basisSet.getBasis(atom);
            val orbitals  = atomBasis.getOrbitals();
            val atombfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} 
                          = new ArrayList[ContractedGaussian{self.at(this)}]();

            for(var orbno:Int=0; orbno<orbitals.size(); orbno++) { 
               val orb = orbitals.get(orbno);
               val typ = orb.getType();
               val pList = plInst.getPowers(typ);

               val coeff:ArrayList[Double]{self.at(this)} = orb.getCoefficients();
               val exps:ArrayList[Double]{self.at(this)}  = orb.getExponents();

               for(var l:Int=0; l<pList.length; l++) {
                  val center = atom.centre;
                  val power = pList(l);
                  val cg:ContractedGaussian{self.at(this)} = new ContractedGaussian(center, power);
                  cg.setIndex(indx++);
                  cg.setIntIndex(intIndx);
              
                  for(var i:Int=0; i<coeff.size(); i++) {
                     cg.addPrimitive(exps.get(i), coeff.get(i));
                  } // end for

                  cg.normalize();

                  basisFunctions.add(cg);
               } // end for

               val am = orb.getAngularMomentum();
               val acg:ContractedGaussian{self.at(this)} = new ContractedGaussian(atom.centre, 
                                                                                  Power(am, 0, 0));
               acg.setIntIndex(intIndx);

               for(var i:Int=0; i<coeff.size(); i++) {
                   acg.addPrimitive(exps.get(i), coeff.get(i));
               } // end for

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
                for(var j:Int=0; j<prms.size(); j++) {
                    val prmj = prms.get(j); 
                    val fact3 = Math.pow(prmj.getExponent(), (2.0*lmn+3.0)/4.0);
                    prmj.setCoefficient(prmj.getCoefficient()*fact1*fact2*fact3);
                } // end for

                val pi3O2By2Tlmn = Math.pow(Math.PI, 1.5) / Math.pow(2.0, lmn); 
                val factC = denom * pi3O2By2Tlmn;
                var factA:Double = 0.0, factB:Double = 0.0;
                var coefExpoSum:Double = 0.0;

                for(var k:Int=0; k<prms.size(); k++) {
                   val prmk = prms.get(k);
                   for(var l:Int=0; l<prms.size(); l++) {
                       val prml = prms.get(l);

                       factA = prmk.getCoefficient() * prml.getCoefficient();
                       factB = Math.pow(prmk.getExponent() + prml.getExponent(), lmn+1.5);
                       coefExpoSum += factA/factB;
                   } // end for
                } // end for 

                val norm = Math.pow(coefExpoSum*factC, -0.5);
                for(var j:Int=0; j<prms.size(); j++) {
                    val prmj = prms.get(j); 
                    prmj.setCoefficient(prmj.getCoefficient()*norm);
                } // end for 
            } // end for
        } // end for
    }

    private def initShellList() {
        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom = molecule.getAtom(atmno);
            val bfs  = atom.getBasisFunctions();
            val nbf  = bfs.size();

            for(var i:Int=0; i<nbf; i++) {
                shellList.addShellPrimitive(bfs.get(i));
            } // end for
        } // end for
        shellList.initPowerList();
    }

    public def getBasisName() = this.basisName;
    public def getBasisFunctions() = basisFunctions;
    public def getShellList() = shellList;
}

