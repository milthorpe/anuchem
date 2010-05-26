/**
 * ConnetivityBuilder.x10
 *
 */
package au.edu.anu.chem;

/**
 * Simple connectivity builder for Molecule object
 * 
 * @author V. Ganesh
 */
public class ConnectivityBuilder[T]{T <: Atom} {
   public def this() {
   }

   /**
    * Build connectivity for the given molecule object
    */
   public def buildConnectivity(mol:Molecule[T]!) {
       val noOfAtoms = mol.getNumberOfAtoms();
       val ai = AtomInfo.getInstance();

       finish foreach(atom in mol.getAtoms()) {
           val conn = new ConnectivitySupport();
           val idx  = atom.getIndex();

           for(var i:Int=0; i<idx; i++) {
              // check for bond between atom and mol.getAtom(i)
              val atomI = mol.getAtom(i);              
              
              if (conn.canFormBond(atom, atomI)) {
                 conn.covalentRadiusSum = ai.getCovalentRadius(atom) + ai.getCovalentRadius(atomI);
                 conn.vdWRadiusSum      = ai.getVdwRadius(atom) + ai.getVdwRadius(atomI);

                 // classify  
                 if (conn.isSingleBondPresent()) {
                    if (conn.isDoubleBondPresent()) {
                       atom.addBond(BondType.DOUBLE_BOND, atomI);
                    } else {
                       atom.addBond(BondType.SINGLE_BOND, atomI);
                    } // end if
                 } // end if
              } // end if 
           }           
       }
   }

   /**
    * detect weak bonds in the molecule object 
    */
   public def detectWeakBonds(mol:Molecule[T]) {
       // TODO:
   }

   /**
    * detect rings (with classification of planar / non-planar) 
    * in the molecule objects
    */
   public def identifyRings(mol:Molecule[T]) {
       // TODO:
   }
}

/** The support class */
class ConnectivitySupport {
    global val BOND_RADIUS = 7.56; // a.u. ~= 4 angstroms
    global val COVALENT_BOND_TOLERANCE = 0.7558903950887472; // a.u. ~= 0.4 angstroms
    global val DOUBLE_BOND_OVERLAP_PERCENTAGE = 0.92;  // 92%
    global val TRIPLE_BOND_OVERLAP_PERCENTAGE = 0.72;  // 72%
   
    var x:Double, y:Double, z:Double;
    var distance:Double;
    var covalentRadiusSum:Double;
    var vdWRadiusSum:Double;

    public def this() { }

    /** Check if the given two atoms can at all form 
        a bond. If yes, selectively store info that will
        be required to later determine the type of bond
      */
    public def canFormBond(a:Atom, b:Atom) : Boolean {
        x = Math.abs(a.centre.i - b.centre.i);
        if (x > BOND_RADIUS) return false;

        y = Math.abs(a.centre.j - b.centre.j);
        if (y > BOND_RADIUS) return false;

        z = Math.abs(a.centre.k - b.centre.k);
        if (z > BOND_RADIUS) return false;

        distance = Math.sqrt(x*x+y*y+z*z);

        return true;
    } 

    /**
     * check for presence of single bonds  
     */
    public def isSingleBondPresent() : Boolean {
        return (((covalentRadiusSum - COVALENT_BOND_TOLERANCE) < distance) 
                && (distance < (covalentRadiusSum + COVALENT_BOND_TOLERANCE)));
    } 

    /**
     * check for presence of double bonds  
     */
    public def isDoubleBondPresent() : Boolean {
        return (distance < (DOUBLE_BOND_OVERLAP_PERCENTAGE * covalentRadiusSum));                
    } // end of method isWeekBondPresent()
}

