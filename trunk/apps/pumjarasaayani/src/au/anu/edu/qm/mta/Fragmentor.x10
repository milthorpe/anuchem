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

import x10.util.ArrayList;
import x10.util.ValRailBuilder;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.BondType;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.ConnectivityBuilder;

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
       for(var i:Int=0; i<noOfAtoms; i++) { 
          mol.getAtom(i).setIndex(i);
          sortedAtomIndices(i) = i;
       } // end for

       // next build the connectivity for this molecule
       val conn = new ConnectivityBuilder[QMAtom]();

       conn.buildConnectivity(mol);
       // conn.detectWeakBonds(mol);
       conn.identifyRings(mol);

       Console.OUT.println("Detected rings: ");
       for(ring in mol.getRings()) {
          Console.OUT.println(ring.isPlanar());
          Console.OUT.println(ring);
       } // end for

       // next start fragmentation, to generate the main fragments

       // step1: atom centered fragments
       generateAtomCenteredFragments(mol, fragList); 
      
       // step2: merge along connectivity path
       mergeAlongConnectivity(mol, sortedAtomIndices, fragList);

       // step3: general merge
       mergeCommonFragments(fragList);

       // step4: purge or expand depending on any rules being broken when a bond is cut
       //        remove dangling atoms, expand double bonds or planar rings 
       finish foreach(fragment in fragList) {
          removeDanglingAtoms(fragment as Fragment!);
       } // finish 

       while(includeMissedAtoms(fragList));

       // step5: generate cardinality fragments 
       new CardinalityExpression().addCardinalityFragments(fragList);

       // step6: add dummy hydrogens, for bonds that are cut
       Console.OUT.println("Adding dummy atoms ...");
       finish foreach(fragment in fragList) {
          addDummyAtoms(fragment as Fragment!);
       } // finish

       // step7: print out general statistics 
       Console.OUT.println("Number of final fragments: " + fragList.size());
       var idx:Int = 0;
       for(frag in fragList) {
           Console.OUT.println("Fragment # " + idx + " : " + frag.getNumberOfTrueAtoms() + ", " 
                               + frag.getNumberOfAtoms() + " [" + frag.cardinalitySign() + "]");
           Console.OUT.println(frag);
           idx++;
       } // end for

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

               if (atom2.getIndex() == atom1.getIndex()) continue;
             
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
 
       Console.OUT.println("Number of bonds for atom " + v + " is " + bonds.size());

       for(var i:Int=0; i<bonds.size(); i++) {
          Console.OUT.println("Processing " + i);
          val bondedAtom = sortedAtomIndices((bonds.get(i).second as QMAtom).getIndex());

          mergeFragmentsCenteredOn(v, bondedAtom, fragList);
          visited(bondedAtom) = true;
          Console.OUT.println("Done processing " + i);
       } //  end for
   }

   /** merge fragments centerd on the given atom indices */
   def mergeFragmentsCenteredOn(a:Int, b:Int, fragList:ArrayList[Fragment]!) {
       Console.OUT.println("Merging " + a + " and " + b + " ...");
       val f1 = findFragmentCenteredOn(a, fragList);
       val f2 = findFragmentCenteredOn(b, fragList);

       if (f1 == null || f2 == null) return;

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
           outer_loop: for(var i:Int=0; i<fragList.size(); i++) {
              val fi = fragList.get(i) as Fragment!;
              for(var j:Int=0; j<fragList.size(); j++) {
                  val fj = fragList.get(j) as Fragment!;
                  doneMerging = (i == fragList.size()-1) && (j == fragList.size()-1);

                  if (i == j) continue;

                  val fifj = fi.union(fj);
                  
                  // check if size constraint is violated, if skip
                  if (fifj.getNumberOfAtoms() > maxFragSize) continue;

                  fragList.remove(fi);
                  fragList.remove(fj);
                  fragList.add(fifj);
                  break outer_loop;
              } // end for
           } // end for            
       } // end while       
   }

   /** include any missed atoms, that violate conditions like double bond breaking
       or breaking rings at inappropriate positions */
   def includeMissedAtoms(fragList:ArrayList[Fragment]!) : Boolean {
       val includedMissedAtoms = Rail.make[Boolean](1, (Int)=>false);

       // TODO: 

       // first see if we are missing any pendant atoms 
       finish foreach(fragment in fragList) {
           val missedAtoms = new ArrayList[QMAtom{self.at(fragment)}]();
           for(atom in fragment.getAtoms()) {
              for(bond in atom.getBonds()) {
                  val bondedAtom = bond.second as QMAtom!;
                  val nBonds = bondedAtom.getBonds().size();

                  if ((nBonds == 1) && (bond.first != BondType.WEAK_BOND)) { 
                      if (!fragment.contains(bondedAtom)) {
                          missedAtoms.add(bondedAtom as QMAtom{self.at(fragment)});
                      } // end if
                  } // end if    
              } // end for
           } // end for

           if (missedAtoms.size() > 0) {
               includedMissedAtoms(0) = true;

               for(matm in missedAtoms) {
                  fragment.addAtom(matm as QMAtom{self.at(fragment)});
               } // end for
           } // end if
       } // end finish

       return includedMissedAtoms(0); 
   }

   /** remove any dangling atoms from a fragment */
   def removeDanglingAtoms(fragment:Fragment!) {
       val ai = AtomInfo.getInstance();
       val atoms = fragment.getAtoms();

       for(atom in atoms) { 
           if (ai.getAtomicNumber(atom) > 2) continue;

           // find connected atoms that are present in the fragment
           val bonds = atom.getBonds(); 
           var nBonds:Int = 0;
           for(bond in bonds) { 
              val bondedAtom = bond.second as QMAtom;

              if (fragment.contains(bondedAtom)) nBonds++;
           } // end for

           if (nBonds == 0) {
               // this is dangling atom
               atoms.remove(atom);
           } // end if 
       } // end for
   }

   /** add dummy atoms to a fragment */
   def addDummyAtoms(fragment:Fragment!) {
       val boundaryAtoms = new ValRailBuilder[QMAtom!]();
 
       val ai = AtomInfo.getInstance(); 
       val hRadius = ai.getCovalentRadius(new QMAtom("H", Point3d(0,0,0), true));

       // find out boundary atoms in this fragment
       for(atom in fragment.getAtoms()) {
           val nBonds  = atom.getBonds().size();
           val nfBonds = fragment.getBondOrder(atom);

           if (nBonds != nfBonds) boundaryAtoms.add(atom);
       } // end for
        
       // then add dummy atoms at appropriate positions
       for(atom in boundaryAtoms.result()) {
           val bonds = atom.getBonds();
           val bondDistance = hRadius + ai.getCovalentRadius(atom as QMAtom!);

           for(bond in bonds) {
              val bondedAtom = bond.second as QMAtom;

              if (!fragment.contains(bondedAtom)) {
                 val vec = (bondedAtom.centre - atom.centre).normalize();
                 val newCenter = Point3d(vec.i*bondDistance + atom.centre.i,
                                         vec.j*bondDistance + atom.centre.j, 
                                         vec.k*bondDistance + atom.centre.k);

                 fragment.addAtom(new QMAtom("H", newCenter, true));
              } // end if
           } // end for
       } // end for
   }

}

