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

import x10x.vector.Point3d;

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
       generateAtomCenteredFragments(mol, fragList); 
      
       // step2: merge along connectivity path
       mergeAlongConnectivity(mol, fragList);

       // step3: general merge

       // step4: purge or expand depending on any rules being broken when a bond is cut

       // step5: add dummy hydrogens, for bonds that are cut

       // step6: print out general statics

       return fragList;
   }

   def generateAtomCenteredFragments(mol:Molecule[QMAtom]!, fragList:ArrayList[Fragment]!) {
       val noOfAtoms = mol.getNumberOfAtoms();

       finish foreach(atom1 in mol.getAtoms()) {
           val aFragment = new Fragment() as Fragment!;
           aFragment.addAtom(atom1);

           for(var i:Int=0; i<noOfAtoms; i++) {
               val atom2 = mol.getAtom(i);
             
               val dist = atom2.centre.distance(atom1.centre);               

               if (dist <= rGoodness) {
                  // include this atom in this fragment
                  aFragment.addAtom(atom2);
               } // end if
           } // end for

           async atomic fragList.add(aFragment); 
       } // finish foreach
   }

   def mergeAlongConnectivity(mol:Molecule[QMAtom]!, fragList:ArrayList[Fragment]!) {
       val noOfAtoms = mol.getNumberOfAtoms();

       for(atom1 in mol.getAtoms()) {
           val bonds = atom1.getBonds();
       
           // use breadth first traversal
       } // finish foreach
   }
}

