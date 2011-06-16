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
package au.edu.anu.chem;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

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
       val ai = AtomInfo.getInstance();

       finish for(atom in mol.getAtoms()) async {
           val conn = new ConnectivitySupport();
           val idx  = atom.getIndex();

           for(var i:Int=0; i<idx; i++) {
              // check for bond between atom and mol.getAtom(i)
              val atomI = mol.getAtom(i);              
              
              if (conn.canFormBond(atom, atomI)) {
                 conn.covalentRadiusSum = ai.getCovalentRadius(atom) + ai.getCovalentRadius(atomI);
                 conn.vdwRadiusSum      = ai.getVdwRadius(atom) + ai.getVdwRadius(atomI);

                 // classify  
                 if (conn.isSingleBondPresent()) {
                    if (conn.isDoubleBondPresent()) {
                       atom.addBond(BondType.DOUBLE_BOND, atomI);
                    } else {
                       atom.addBond(BondType.SINGLE_BOND, atomI);
                    } // end if
                 } // end if
              } // end if 
           } // end for           
       } // finish
   }

   /**
    * detect weak bonds in the molecule object 
    */
   public def detectWeakBonds(mol:Molecule[T]) {
       // TODO:
       val noOfAtoms = mol.getNumberOfAtoms();
       val ai = AtomInfo.getInstance();
       val WEAK_BOND_AXIS_REJECTION_FACTOR = 0.1;
       val WEAK_BOND_ANGLE_TOLERANCE = 1.222;

       finish for(atom in mol.getAtoms()) async {
           val conn = new ConnectivitySupport();
           val idx  = atom.getIndex();

           for(var i:Int=0; i<idx; i++) {
              // check for bond between atom and mol.getAtom(i)
              val atomI = mol.getAtom(i);

              if (conn.canFormBond(atom, atomI)) {
                 conn.covalentRadiusSum = ai.getCovalentRadius(atom) + ai.getCovalentRadius(atomI);
                 conn.vdwRadiusSum      = ai.getVdwRadius(atom) + ai.getVdwRadius(atomI);

                 if (conn.isWeekBondPresent() && (!conn.isSingleBondPresent())) {
                    // TODO: make it specific to type of atoms
                    val bonds1 = atom.getBonds();
                    val bonds2 = atomI.getBonds();

                    if ((bonds1.size()==0) || (bonds2.size()==0)) {
                       atom.addBond(BondType.WEAK_BOND, atomI);
                       continue;
                    } // end if

                    val axis1 = computeAxis(atom); 
                    val axis2 = computeAxis(atomI); 

                    if ((axis1==Vector3d.NULL) || (axis2==Vector3d.NULL)) continue;
       
                    if (axis1.magnitude() < WEAK_BOND_AXIS_REJECTION_FACTOR
                        || axis2.magnitude() < WEAK_BOND_AXIS_REJECTION_FACTOR) {
                        continue;
                    } // end if

                    val v12 = atomI.centre - atom.centre;
                    val v21 = atom.centre  - atomI.centre;
                    
                    val angle1 = axis1.angleWith(v12);
                    val angle2 = axis2.angleWith(v21);

                    if ((angle1 < WEAK_BOND_ANGLE_TOLERANCE) 
                        && (angle2 < WEAK_BOND_ANGLE_TOLERANCE)) {
                       atom.addBond(BondType.WEAK_BOND, atomI);
                    } // end if
                 } // end if
              } // end if
           } // end for
       } // finish
   }

   private def computeAxis(atom:Atom) : Vector3d {
       val bonds = atom.getBonds();
       var axis:Vector3d  = Vector3d.NULL;

       for(bond in bonds) {
          if ((bond.first == BondType.SINGLE_BOND) || (bond.first == BondType.DOUBLE_BOND)) {
             axis = axis + (atom.centre - bond.second.centre);
          } // end if
       } // end for

       return axis;
   }

   private val WHITE = 0; 
   private val GRAY = 1; 
   private val BLACK = 2;    
   private val TORSSIAN_ANGLE_TOLERANCE = 0.08726646259971647; // radians

   /**
    * detect rings (with classification of planar / non-planar) 
    * in the molecule objects
    */
   public def identifyRings(mol:Molecule[T]) {
       val noOfAtoms = mol.getNumberOfAtoms();
       
       val color  = new Array[Int](noOfAtoms, (Int)=>WHITE);
       val parent = new Array[Int](noOfAtoms, (Int)=>-1);

       // detect rings
       for(var i:Int=0; i<noOfAtoms; i++) {
          if (color(i) == WHITE) {
             traverseAndRecordRing(mol, i, color, parent);   
          } // end if 
       } // end for
       
       // TODO: remove subsets

       // set planarity
       for(ring in mol.getRings()) {
          val n = ring.getSize();

          if (n < 4) {
             ring.setPlanar(true);
             continue;
          } // end if

          val points = new Array[Point3d](n);
          
          for(var i:Int=0; i<n; i++) {
             points(i) = ring.getAtom(i).centre; 
          } // end for

          for(var i:Int=0; i<n-1; i++) {
             val a12 = points(i+1) - points(i);
             val i2 = ((i+2 >= n) ? 0 : (i + 2));
           
             val a23 = points(i2) - points(i+1);            
             val i3 = ((i+3 >= n) ? 0 : (i + 3));
            
             val a34 = points(i3) - points(i2);
            
             val n1 = a12.cross(a23);
             val n2 = a23.cross(a34);

             if (n1.angleWith(n2) > TORSSIAN_ANGLE_TOLERANCE) {
                ring.setPlanar(false); // not planar
                break;
             } // end if 
            
          } // end for
       } // end for
   }

   private def traverseAndRecordRing(mol:Molecule[T], v:Int, color:Rail[Int], parent:Rail[Int]) {
        color(v) = GRAY;  // this vertex is to be processed now
        
        val bonds = mol.getAtom(v).getBonds();
        var atomIndex:Int;
        
        // traverse the strongly connected ones only
        for(bond in bonds) {
            atomIndex = bond.second.getIndex();
            
            if (bond.first != BondType.WEAK_BOND || bond.first != BondType.NO_BOND) {
                if (color(atomIndex) == WHITE) {   // not yet traversed?
                    parent(atomIndex) = v;
                    traverseAndRecordRing(mol, atomIndex, color, parent); // recursive call
                } else if (color(atomIndex) == BLACK) { // found a ring!
                    val theRing = new Ring[T]();
                    
                    theRing.addAtom(mol.getAtom(atomIndex));
                    
                    var vertex:Int = atomIndex;
                    
                    // record the cycle
                    while(true) {
                        vertex = parent(vertex);
                        
                        if (vertex == -1) { 
                            break;
                        } else if (vertex == v) {
                            theRing.addAtom(mol.getAtom(vertex));
                            break;
                        } else {
                            theRing.addAtom(mol.getAtom(vertex));
                        } // end if
                    } // end while  
                    
                    // record the ring
                    mol.addRing(theRing);
                } // end if
            } // end if
        } // end while
        
        // this node has been completely traversed, if we hit this node again
        // we have a cycle!
        color(v) = BLACK;
   }

    /** The support class */
    static class ConnectivitySupport {
        static val BOND_RADIUS = 7.56; // a.u. ~= 4 angstroms
        static val COVALENT_BOND_TOLERANCE = 0.7558903950887472; // a.u. ~= 0.4 angstroms
        static val WEAK_BOND_TOLERANCE_LOWER = 0.1889725987721868; // a.u. ~= 0.1 angstroms
        static val WEAK_BOND_TOLERANCE_UPPER = 1.9842122871079613; // a.u. ~= 1.05 angstroms
        static val DOUBLE_BOND_OVERLAP_PERCENTAGE = 0.92;  // 92%
        static val TRIPLE_BOND_OVERLAP_PERCENTAGE = 0.72;  // 72%
       
        var x:Double, y:Double, z:Double;
        var distance:Double;
        var covalentRadiusSum:Double;
        var vdwRadiusSum:Double;

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

        /**
         * method to check the presence of weak bond, using distance criterion
         */
        public def isWeekBondPresent() : Boolean {
            return ((distance < (vdwRadiusSum - WEAK_BOND_TOLERANCE_LOWER)
                 && (vdwRadiusSum - WEAK_BOND_TOLERANCE_UPPER) < distance));
        } 
    }
}


