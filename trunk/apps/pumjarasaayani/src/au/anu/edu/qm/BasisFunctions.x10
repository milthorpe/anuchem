/**
 * BasisFunctions.x10
 *
 * Generates basis functions for a given Molecule object 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.ArrayList;
import au.edu.anu.chem.Molecule;

public class BasisFunctions { 
    global var molecule:Molecule[QMAtom]{self.at(this)};
    global var basisName:String;
    global var basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    global var shellList:ShellList{self.at(this)};

    public def this() { }

    public def make(mol:Molecule[QMAtom]{self.at(this)}, basNam:String, basisDir:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian{self.at(this)}](); 
        initBasisFunctions(basisDir);

        shellList = new ShellList();
        shellList.make();
        initShellList();
    } 

    private def initBasisFunctions(basisDir:String) {
        var basisSet:BasisSet{self.at(this)} = new BasisSet(); 
        basisSet.make(basisName, basisDir);
        var indx:Int = 0;
        var intIndx:Int = 0;

        for(var atmno:Int=0; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom      = molecule.getAtom(atmno);
            val atomBasis = basisSet.getBasis(atom);
            val orbitals  = atomBasis.getOrbitals();
            val atombfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} 
                          = new ArrayList[ContractedGaussian{self.at(this)}]();

            for(var orbno:Int=0; orbno<orbitals.size(); orbno++) { 
               val orb = orbitals.get(orbno);
               val typ = orb.getType();
               val pList = PowerList.getInstance().getPowers(typ);
               val pListSiz = pList.region.max(0); 

               for(var l:Int=0; l<pListSiz; l++) {
                  var cg:ContractedGaussian{self.at(this)} = new ContractedGaussian();
                  cg.make(atom.centre, pList(l));
                  cg.setIndex(indx++);
                  cg.setIntIndex(intIndx);
              
                  val coeff:ArrayList[Double]{self.at(this)} = orb.getCoefficients();
                  val exps:ArrayList[Double]{self.at(this)}  = orb.getExponents();

                  for(var i:Int=0; i<coeff.size(); i++) {
                     cg.addPrimitive(exps.get(i), coeff.get(i));
                     intIndx++;
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
    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;
    public def getShellList() : ShellList{self.at(this)} = shellList;
}

