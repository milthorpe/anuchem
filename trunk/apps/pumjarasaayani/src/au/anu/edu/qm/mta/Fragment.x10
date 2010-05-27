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

     private var centeredOn:Int;

     private var energy:Double;

     private var cardinalitySign:Int;
     
     public def this() {
          centeredOn = -1;
          cardinalitySign = 1;
     }

     public def centeredOn() = centeredOn;
     public def centeredOn(atmIdx:Int) { centeredOn = atmIdx; }

     public def energy(ene:Double) { energy = ene; }
     public def energy() = energy;

     public def cardinalitySign(sign:Int) { cardinalitySign = sign; }
     public def cardinalitySign() = cardinalitySign;

     public def intersection(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

          var foundAtom:Boolean;

          for(atom1 in getAtoms()) {
             val idx = atom1.getIndex();
             foundAtom = false;
             for(atom2 in frag.getAtoms()) {
                 if (atom2.getIndex() == idx) {
                    foundAtom = true; break;
                 } // end if
             } // end for

             if (foundAtom) newFrag.addAtom(atom1);
          } // end for

          return newFrag;
     }

     public def union(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

          var foundAtom:Boolean;

          for(atom1 in frag.getAtoms()) { 
             newFrag.addAtom(atom1);
          } // end for

          for(atom1 in getAtoms()) {
             val idx = atom1.getIndex();
             foundAtom = false;
             for(atom2 in newFrag.getAtoms()) {
                 if (atom2.getIndex() == idx) {
                    foundAtom = true; break;
                 } // end if
             } // end for

             if (!foundAtom) newFrag.addAtom(atom1);
          } // end for

          return newFrag;
     } 

     public def contains(atm:QMAtom) : Boolean {
          val idx = atm.getIndex();

          for(atom in getAtoms()) {
             if ((atom as QMAtom).getIndex() == idx) return true;
          } // end for

          return false;
     }

     public def getBondOrder(atm:QMAtom) : Int {
          if (!contains(atm)) return 0;

          val bonds = atm.getBonds();
          var nBonds:Int = 0;
          for(bond in bonds) {
             if (contains(bond.second as QMAtom)) nBonds++;
          } // end for

          return nBonds;
     }
}

