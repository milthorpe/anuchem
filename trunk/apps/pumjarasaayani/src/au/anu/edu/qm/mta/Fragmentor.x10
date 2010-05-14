/**
 * Fragmentor.x10
 *
 * Basic implementation of MTA in X10, for large scale HF calculation. 
 * Lots of code derived from MeTA Studio and MTA-GAMESS written by me ;)
 * This is a straight implementation, so no plugins are provided as 
 * with MeTA Studio.
 *
 * TODO: Currently this is a very basic implementation, does not take 
 * care of all the conditions as implemented in MTA-GAMESS.
 *
 * Primary reference:
 *   [MTA06] V. Ganesh, R. K. Dongare, P. Balanarayan, and S. R. Gadre, 
 *           J. Chem. Phys. 125, 104109, 2006.
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm.mta;

import x10.util.ArrayList;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.ConnectivityBuilder;

public class Fragmentor {

   val rGoodness:Double;  // note that this is in a.u. as opposed to the original code!
   val maxFragSize:Int;

   public def this(rGoodness:Double, maxFragSize:Int) {
       this.rGoodness = rGoodness;
       this.maxFragSize = maxFragSize;
   }

   public def fragment(mol:Molecule[QMAtom]!) : ArrayList[Fragment]! {
       val fragList = new ArrayList[Fragment]() as ArrayList[Fragment]!;

       // TODO: reorder the atom indices

       // next build the connectivity for this molecule
       val conn = new ConnectivityBuilder[QMAtom]();

       conn.buildConnectivity(mol);
       conn.detectWeakBonds(mol);
       conn.identifyRings(mol);

       // next start fragmentation, to generate the main fragments

       // step1: atom centered fragments
       
      
       // step2: merge along connectivity path

       // step3: general merge

       // step4: add dummy hydrogens, for bonds that are cut

       // step5: print out general statics

       return fragList;
   }

   def generateAtomCenteredFragments(mol:Molecule[QMAtom]!, fragList:ArrayList[Fragment]!) {
       val noOfAtoms = mol.getNumberOfAtoms();

       finish foreach(atom1 in mol.getAtoms()) {
           val aFragment = new Fragment() as Fragment!;
           aFragment.addAtom(atom1);

           for(var i:Int=0; i<noOfAtoms; i++) {
               val atom2 = mol.getAtom(i);
             
               val dist = atom2.center.distance(atom1);               
           } // end for

           async atomic fragList.add(aFragment); 
       } // finish foreach
   }
}

