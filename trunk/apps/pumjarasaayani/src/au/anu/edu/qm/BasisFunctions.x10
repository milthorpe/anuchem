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
    global val molecule:Molecule[QMAtom]{self.at(this)};
    global val basisName:String;
    global val basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    global val shellList:ShellList{self.at(this)};

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

               val coeff:ArrayList[Double]{self.at(this)} = orb.getCoefficients();
               val exps:ArrayList[Double]{self.at(this)}  = orb.getExponents();

               val norm = Rail.make[Double](coeff.size(), (Int)=>1.0);

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
                                                                                  new Power(am, 0, 0));
               acg.setIntIndex(intIndx);

               for(var i:Int=0; i<coeff.size(); i++) {
                   acg.addPrimitive(exps.get(i), coeff.get(i));
               } // end for

               acg.normalize();

               for(var i:Int=0; i<coeff.size(); i++) {
                   val pg = acg.getPrimitive(i);
                   pg.setCoefficient(coeff.get(i) * acg.getNormalization() * pg.getNormalization());
               } // end for

               atombfs.add(acg);
               intIndx += ((am+1)*(am+2)/2);
            } // end for          

            atom.setBasisFunctions(atombfs);  
        } // end for
    }

    private def initShellList() {
        for(cg in basisFunctions) shellList.addShellPrimitive(cg);
        shellList.initPowerList();
    }

    public def getBasisName() = this.basisName;
    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;
    public def getShellList() : ShellList{self.at(this)} = shellList;
}

