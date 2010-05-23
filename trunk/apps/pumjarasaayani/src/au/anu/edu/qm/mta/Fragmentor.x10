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

   /** creates a new instance of Fragmentor.
       See [MTA06] for explanation on rGoodness and maxFragSize */
   public def this(rGoodness:Double, maxFragSize:Int) {
       this.rGoodness = rGoodness;
       this.maxFragSize = maxFragSize;
   }

   /** fragment the molecule into overlapping fragments */
   public def fragment(mol:Molecule[QMAtom]!) : ArrayList[Fragment]! {
       val fragList = new ArrayList[Fragment]() as ArrayList[Fragment]!;

       val noOfAtoms = mol.getNumberOfAtoms();

       // TODO: reorder the atom indices, depending on nearness to the 
       // center of mass of the molecule so as to ensure a more 
       // consistant ordering that is indipendent of the input order
       val sortedAtomIndices = Rail.make[Int](noOfAtoms) as Rail[Int]!;
       for(var i:Int=0; i<noOfAtoms; i++) sortedAtomIndices(i) = i;

       // next build the connectivity for this molecule
       val conn = new ConnectivityBuilder[QMAtom]();

       conn.buildConnectivity(mol);
       conn.detectWeakBonds(mol);
       conn.identifyRings(mol);

       // next start fragmentation, to generate the main fragments

       // step1: atom centered fragments
       generateAtomCenteredFragments(mol, fragList); 
      
       // step2: merge along connectivity path
       mergeAlongConnectivity(mol, sortedAtomIndices, fragList);

       // step3: general merge
       mergeCommonFragments(mol);

       // step4: purge or expand depending on any rules being broken when a bond is cut
       //        remove dangling atoms, expand double bonds or planar rings 

       // step5: add dummy hydrogens, for bonds that are cut

       // step6: print out general statistics 

       return fragList;
   }

   /** generate atom centered fragments with minimum r-goodness as requeated */
   def generateAtomCenteredFragments(mol:Molecule[QMAtom]!, fragList:ArrayList[Fragment]!) {
       val noOfAtoms = mol.getNumberOfAtoms();

       finish foreach(atom1 in mol.getAtoms()) {
           val aFragment = new Fragment() as Fragment!;
           aFragment.centeredOn(atom1.getIndex());
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

   /** merge fragments along the connectivity path */
   def mergeAlongConnectivity(mol:Molecule[QMAtom]!, sortedAtomIndices:Rail[Int]!, fragList:ArrayList[Fragment]!) {
       val noOfAtoms = mol.getNumberOfAtoms(); 
       val visited = Rail.make[Boolean](noOfAtoms, (Int)=>false) as Rail[Boolean]!;        
        
       for(var i:Int=0; i<noOfAtoms; i++) {
           val idx = sortedAtomIndices(i);
           val atom1 = mol.getAtom(idx);
           val bonds = atom1.getBonds();

           if (!visited(idx)) {
              visited(idx) = true;
   
              traverseAndMergeFragments(idx, sortedAtomIndices, mol, visited, fragList); 
           } // end if
       } // finish foreach
   }

   /** simple traversal and merge */
   def traverseAndMergeFragments(v:Int, sortedAtomIndices:Rail[Int]!, mol:Molecule[QMAtom]!, visited:Rail[Boolean]!, fragList:ArrayList[Fragment]!) {
       val bonds = mol.getAtom(v).getBonds();
 
       for(var i:Int=0; i<bonds.size(); i++) {
          val bondedAtom = sortedAtomIndices((bonds.get(i).second as QMAtom).getIndex());

          mergeFragmentsCenteredOn(v, bondedAtom, fragList);
          visited(bondedAtom) = true;
       } //  end for
   }

   /** merge fragments centerd on the given atom indices */
   def mergeFragmentsCenteredOn(a:Int, b:Int, fragList:ArrayList[Fragment]!) {
       val f1 = findFragmentCenteredOn(a, fragList);
       val f2 = findFragmentCenteredOn(b, fragList);

       val f1f2 = f1.union(f2);
       // check if size constraint is violated, if so exit
       if (f1f2.getNumberOfAtoms() > maxFragSize) return;

       // else remove f1 and f2 from fragList and add f1f2 instead
       fragList.remove(f1);
       fragList.remove(f2);

       f1f2.centeredOn(f1.centeredOn());
       fragList.add(f1f2);
   } 

   /** return the atom centered fragment */
   def findFragmentCenteredOn(idx:Int, fragList:ArrayList[Fragment]!) : Fragment! {
       for(frag in fragList) {
           if (frag.centeredOn() == idx) return (frag as Fragment!);
       } // end for

       return null;  // should never come here!
   }

   /** general merge procedure, keep merging until the maxFragmentSize criteria can not be met */
   def mergeCommonFragments(fragList:ArrayList[Fragment]!) {
       var doneMerging:Boolean = false;

       while(!doneMerging) {
           val noOfFragments = fragList.size();

           for(var i:Int=0; i<noOfFragments; i++) {
              for(var j:Int=0; j<noOfFragments; j++) {
                  doneMerging = (i == noOfFragments) && (j == noOfFragments);

                  if (i == j) continue;

                  // TODO:
              } // end for
           } // end for            
       } // end while       
   }

}

