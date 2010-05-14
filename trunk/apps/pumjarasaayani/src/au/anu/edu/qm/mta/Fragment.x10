/**
 * Fragment.x10
 * 
 * Represents a single fragment.
 *
 * @author: V.Ganesh
 */


package au.anu.edu.qm.mta;

import x10.util.ArrayList;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.Molecule;

public class Fragment extends Molecule[QMAtom] {

     val dummyAtoms:ArrayList[QMAtom];
     
     public def this() {
         dummyAtoms = new ArrayList[QMAtom]();
     }

     public def addDummyAtom(dummyAtom:QMAtom!) {
         dummyAtoms.add(dummyAtom);
     }

     public def intersection(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

          return newFrag;
     }

     public def union(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

          return newFrag;
     } 
}

