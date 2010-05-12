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
   public def buildConnectivity(mol:Molecule[T]) {
       val noOfAtoms = mol.getNumberOfAtoms();
       val ai = new AtomInfo();

       finish foreach(atom in mol.getAtoms()) {
           val conn = new ConnectivitySupport(ai);

           for(var i:Int=0; i<noOfAtoms; i++) {
              // check for bond between atom and mol.getAtom(i)
              val atomI = mol.getAtom(i);
              
              if (conn.canFormBond(atom, atomI)) {
                 // classify  
                 if (conn.isSingleBondPresent()) {
                    if (conn.isDoubleBondPresent()) {
                    } // end if
            
                    atom.addBond(BondType.SINGLE_BOND, atomI);
                    atomI.addBond(BondType.SINGLE_BOND, atom);
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

    val ai:AtomInfo;

    def this(ai:AtomInfo) { this.ai = ai; }

    /** Check if the given two atoms can at all form 
        a bond. If yes, selectively store info that will
        be required to later determine the type of bond
      */
    def canFormBond(a:Atom, b:Atom) : Boolean {
        x = Math.abs(a.centre.i - b.centre.i);
        if (x > BOND_RADIUS) return false;

        y = Math.abs(a.centre.j - b.centre.j);
        if (y > BOND_RADIUS) return false;

        z = Math.abs(a.centre.k - b.centre.k);
        if (z > BOND_RADIUS) return false;

        distance = Math.sqrt(x*x+y*y+z*z);

        covalentRadiusSum = ai.getCovalentRadius(a) + ai.getCovalentRadius(b);
        vdWRadiusSum = ai.getVdwRadius(a) + ai.getVdwRadius(b);

        return true;
    } 

    /**
     * check for presence of single bonds  
     */
    def isSingleBondPresent() : Boolean {
        return (((covalentRadiusSum - COVALENT_BOND_TOLERANCE) < distance) 
                && (distance < (covalentRadiusSum + COVALENT_BOND_TOLERANCE)));
    } 

    /**
     * check for presence of double bonds  
     */
    def isDoubleBondPresent() : Boolean {
        return (distance < (DOUBLE_BOND_OVERLAP_PERCENTAGE * covalentRadiusSum));                
    } // end of method isWeekBondPresent()
}

