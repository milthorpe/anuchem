/**
 * Fragment.x10
 * 
 * Represents a single fragment.
 *
 * @author: V.Ganesh
 */


package au.anu.edu.qm.mta;

import x10.util.Pair;
import x10.util.ArrayList;
import x10.util.ValRailBuilder;

import x10x.vector.Point3d;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.Molecule;

public class Fragment extends Molecule[QMAtom] {

     global val dummyAtoms:ArrayList[QMAtom];
     
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

          var foundAtom:Boolean;

          for(atom1 in frag.getAtoms()) { 
             // TODO: contains pattern, move out
             val idx = atom1.getIndex();

             foundAtom = false;
             for(atom2 in getAtoms()) {
                 if (atom2.getIndex() == idx) {
                    foundAtom = true; break;
                 } // end if
             } // end for

             if (!foundAtom) addAtom(atom1);
          } // end for

          return newFrag;
     } 

     public safe def getCoords() : ValRail[Pair[String,Point3d]] {
          val noOfAtoms = getNumberOfAtoms() + dummyAtoms.size();
          val coords = new ValRailBuilder[Pair[String,Point3d]](noOfAtoms);

          for(atom in getAtoms()) {
              coords.add(Pair[String,Point3d](atom.symbol, atom.centre));
          } // end for
          
          for(atom in dummyAtoms) {
              coords.add(Pair[String,Point3d](atom.symbol, atom.centre));
          } // end for

          return coords.result();
     }
}

