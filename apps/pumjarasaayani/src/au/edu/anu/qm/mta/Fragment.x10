/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.edu.anu.qm.mta;

import x10.util.ArrayList;

import x10x.vector.Point3d;

import au.edu.anu.qm.QMAtom;
import au.edu.anu.chem.Molecule;

/**
 * Fragment.x10
 *
 * Represents a single fragment.
 *
 * @author: V.Ganesh
 */
public class Fragment extends Molecule[QMAtom] {

     public var centreedOn:Int;

     public var energy:Double;

     public var cardinalitySign:Int;
     
     public def this() {
          centreedOn = -1;
          cardinalitySign = 1;
     }

     public def getNumberOfTrueAtoms() : Int {
          var nAtoms:Int = 0;
 
          for(atom in getAtoms()) 
             if (!(atom as QMAtom).isDummy()) nAtoms++;

          return nAtoms;
     }

     public def intersection(frag:Fragment) : Fragment {
          val newFrag = new Fragment();

          var foundAtom:Boolean;

          for(atom1 in getAtoms()) {
             if ((atom1 as QMAtom).isDummy()) continue;

             val idx = atom1.getIndex();
             foundAtom = false;
             for(atom2 in frag.getAtoms()) {
                 if ((atom2 as QMAtom).isDummy()) continue;

                 if (atom2.getIndex() == idx) {
                    foundAtom = true; break;
                 } // end if
             } // end for

             if (foundAtom) newFrag.addAtom(atom1);
          } // end for

          return newFrag;
     }

     public def union(frag:Fragment) : Fragment {
          val newFrag = new Fragment();

          var foundAtom:Boolean;

          for(atom1 in frag.getAtoms()) { 
             if ((atom1 as QMAtom).isDummy()) continue;
             newFrag.addAtom(atom1);
          } // end for

          for(atom1 in getAtoms()) {
             if ((atom1 as QMAtom).isDummy()) continue;
             val idx = atom1.getIndex();
             foundAtom = false;
             for(atom2 in newFrag.getAtoms()) {
                 if ((atom2 as QMAtom).isDummy()) continue;
                 if (atom2.getIndex() == idx) {
                    foundAtom = true; break;
                 } // end if
             } // end for

             if (!foundAtom) newFrag.addAtom(atom1);
          } // end for

          return newFrag;
     } 

     public def equals(frag:Fragment) : Boolean {          

          if (getNumberOfTrueAtoms() != frag.getNumberOfTrueAtoms()) return false;

          for(atom1 in frag.getAtoms()) {
              val atm = atom1 as QMAtom;
              if (atm.isDummy()) continue;
              if (!contains(atm)) return false;
          } // end for

          return true;
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

