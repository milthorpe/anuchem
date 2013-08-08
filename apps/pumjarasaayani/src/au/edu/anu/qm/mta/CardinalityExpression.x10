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

/**
 * CardinalityExpression.x10
 *
 * Generate cardinality expression for patching energy, or other related properties
 * given a fragmentation scheme.
 *
 * Primary reference:
 *   V. Ganesh, R. K. Dongare, P. Balanarayan, and S. R. Gadre,
 *   J. Chem. Phys. 125, 104109, 2006.
 *
 * @author: V.Ganesh
 */
public class CardinalityExpression {

   private var mergedIntoPreviousTerm:Boolean;

   public def this() {
   }

   /** add cardinality fragments */
   public def addCardinalityFragments(fragList:ArrayList[Fragment]) {
       val noOfFragments = fragList.size() as Int;
       val combs = new Rail[Int](noOfFragments, -1n);
       val cfList = new ArrayList[Fragment]();
       var l:Int, pos:Int, m:Int;

       for(var i:Int=0n; i<noOfFragments; i++) {
          for(var j:Int=i+1n; j<noOfFragments; j++) {
             combs(0) = i; combs(1) = j;

             if (!computeIntersections(fragList, cfList, combs, 1n, getSignOfTerm(1n))) continue;

             if (noOfFragments <= 2n || combs(1) == noOfFragments-1n) continue;

             l=1n; pos=1n; m=pos+1n;

             while(true) {
                 for (k in (combs(pos)+1n)..(noOfFragments-1n)) {
                     combs(m) = k;
                        
                     if (m > 0n) {
                         if (combs(m) <= combs(m-1)) break;
                     } // end if

                     if (!computeIntersections(fragList, cfList, combs, m, getSignOfTerm(m))) {
                         if (combs(m) == noOfFragments-1n) break;
                         m--;
                     } // end  if
                        
                     m++;
                 } // end for

                 pos=noOfFragments-l; m=pos;
                    
                 if (combs(pos) == noOfFragments-1n) l++;
                 else                               l=1n;
                    
                 if (pos==1n) break;
             } // end while
          } // end for
       } // end for

       Console.OUT.println("No. of cardinality fragments: " + cfList.size());

       for(cf in cfList) {        
          if (cf.cardinalitySign != 0n) fragList.add(cf);
       } // end for
   }

   private def getSignOfTerm(n:Int) : Int {
       val n1 = n+1n;

       if (n1%2n == 0n) return -1n;
       else             return 1n;
   }

   private def computeIntersections(fragList:ArrayList[Fragment], cfList:ArrayList[Fragment], combs:Rail[Int], nTerms:Int, sign:Int) : Boolean {
       var fInter:Fragment = fragList.get(combs(0));       

       for(var i:Int=1n; i<=nTerms; i++) {
           fInter = fInter.intersection(fragList.get(combs(i)));

           if (fInter.getNumberOfTrueAtoms() == 0L) return false;
       } // end for

       addFragmentToList(fragList, cfList, fInter, sign);
       return true;
   }

   private def addFragmentToList(fragList:ArrayList[Fragment], cfList:ArrayList[Fragment], fragment:Fragment, sign:Int) {
       mergedIntoPreviousTerm = false;

       for(frag in cfList) {
          if (frag.equals(fragment)) {
             mergedIntoPreviousTerm = true;
             frag.cardinalitySign = frag.cardinalitySign + sign;
             break;
          } // end if
       } // end for

       if (!mergedIntoPreviousTerm) { 
           fragment.cardinalitySign = sign;
           cfList.add(fragment);
       } // end if
   }
}

