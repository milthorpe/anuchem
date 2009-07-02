package au.anu.edu.qm;

import x10.util.*;

public class BasisFunctions { 
    var molecule:Molecule;
    var basisName:String;
    var basisFunctions:ArrayList[ContractedGaussian];

    public def this(mol:Molecule, basNam:String) { 
        this.molecule  = mol;
        this.basisName = basNam;

        basisFunctions = new ArrayList[ContractedGaussian](); 
        initBasisFunctions();
    } 

    private def initBasisFunctions() : void {
        val basisSet:BasisSet = new BasisSet(basisName);

        for(atom in molecule.getAtoms()) {
            val orbitals = basisSet.getBasis(atom).getOrbitals();

            for(orb in orbitals) { 
               val pList = PowerList.getInstance().getPowers(orb.getType());
              
               for(var l:Int=0; l<pList.region.max(0); l++) {
                  var cg:ContractedGaussian = new ContractedGaussian(atom, pList(l));   
                  val coeff = orb.getCoefficients();
                  val exps  = orb.getExponents();

                  for(var i:Int=0; i<coeff.size(); i++) {
                     cg.addPrimitive(exps.get(i), coeff.get(i));
                  } // end for

                  cg.normalize();
                  basisFunctions.add(cg);
               } // end for
            } // end for            
        } // end for
    }

    public def getBasisName() : String = this.basisName;

    public def getBasisFunctions() : ArrayList[ContractedGaussian] {    
        return basisFunctions;
    }
}

