/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm.mta;

import x10.util.Pair;
import x10.util.ArrayList;
import x10.util.ValRailBuilder;

import x10x.vector.Point3d;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.Molecule;

/**
 * Fragment.x10
 *
 * Represents a single fragment.
 *
 * @author: V.Ganesh
 */
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

     public def getNumberOfTrueAtoms() : Int {
          var nAtoms:Int = 0;
 
          for(atom in getAtoms()) 
             if (!(atom as QMAtom).isDummy()) nAtoms++;

          return nAtoms;
     }

     public def intersection(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

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

     public def union(frag:Fragment!) : Fragment! {
          val newFrag = new Fragment() as Fragment!;

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

     public def equals(frag:Fragment!) : Boolean {          

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

