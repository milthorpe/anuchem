/**
 * BasisFunctions.x10
 *
 * Generates basis functions for a given Molecule object 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class BasisFunctions { 
    var molecule:Molecule;
    var basisName:String;
    var basisFunctions:ArrayList[ContractedGaussian];
    var shellList:ShellList;

    public def this() { }

    public def make(mol:Molecule, basNam:String, basisDir:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian](); 
        initBasisFunctions(basisDir);

        shellList = new ShellList();
        initShellList();
    } 

    private def initBasisFunctions(basisDir:String) {
        var basisSet:BasisSet = new BasisSet(); 
        basisSet.make(basisName, basisDir);
        var indx:Int = 0;

        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom      = molecule.getAtom(atmno);
            val atomBasis = basisSet.getBasis(atom);
            val orbitals  = atomBasis.getOrbitals();
            val atombfs   = new ArrayList[ContractedGaussian]();

            for(var orbno:Int=0; orbno<orbitals.size(); orbno++) { 
               val orb = orbitals.get(orbno);
               val pList = PowerList.getInstance().getPowers(orb.getType());
               val pListSiz = pList.region.max(0); 

               for(var l:Int=0; l<pListSiz; l++) {
                  var cg:ContractedGaussian = new ContractedGaussian(atom, pList(l));   
                  cg.setIndex(indx++);
              
                  val coeff:ArrayList[Double] = orb.getCoefficients();
                  val exps:ArrayList[Double]  = orb.getExponents();

                  for(var i:Int=0; i<coeff.size(); i++) {
                     cg.addPrimitive(exps.get(i), coeff.get(i));
                  } // end for

                  cg.normalize();
                  atombfs.add(cg);
                  basisFunctions.add(cg);
               } // end for
            } // end for          

            atom.setBasisFunctions(atombfs);  
        } // end for
    }

    private def initShellList() {
        for(cg in basisFunctions) shellList.addShellPrimitive(cg);
    }

    public def getBasisName() = this.basisName;
    public def getBasisFunctions() = basisFunctions;
    public def getShellList() = shellList;
}

