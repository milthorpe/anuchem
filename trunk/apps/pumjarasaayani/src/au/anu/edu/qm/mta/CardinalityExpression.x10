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

package au.anu.edu.qm.mta;

import x10.util.ArrayList;
import x10.util.ValRailBuilder;

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
   public def addCardinalityFragments(fragList:ArrayList[Fragment!]!) {
       val noOfFragments = fragList.size();
       val combs = Rail.make[Int](noOfFragments, (Int)=>-1);
       val cfList = new ArrayList[Fragment!]();
       var l:Int, pos:Int, m:Int;

       for(var i:Int=0; i<noOfFragments; i++) {
          for(var j:Int=i+1; j<noOfFragments; j++) {
             combs(0) = i; combs(1) = j;

             if (!computeIntersections(fragList, cfList, combs, 1, getSignOfTerm(1))) continue;

             if (noOfFragments <= 2 || combs(1) == noOfFragments-1) continue;

             l=1; pos=1; m=pos+1;

             while(true) {
                 for(var k:Int=combs(pos)+1; k<noOfFragments; k++) {
                     combs(m) = k;
                        
                     if (m > 0) {
                         if (combs(m) <= combs(m-1)) break;
                     } // end if

                     if (!computeIntersections(fragList, cfList, combs, m, getSignOfTerm(m))) {
                         if (combs(m) == noOfFragments-1) break;
                         m--;
                     } // end  if
                        
                     m++;
                 } // end for

                 pos=noOfFragments-l; m=pos;
                    
                 if (combs(pos) == noOfFragments-1) l++;
                 else                               l=1;
                    
                 if (pos==1) break;
             } // end while
          } // end for
       } // end for

       Console.OUT.println("No. of cardinality fragments: " + cfList.size());

       for(cf in cfList) {        
          if (cf.cardinalitySign != 0) fragList.add(cf);
       } // end for
   }

   private def getSignOfTerm(n:Int) : Int {
       val n1 = n+1;

       if (n1%2 == 0) return -1;
       else           return 1;
   }

   private def computeIntersections(fragList:ArrayList[Fragment!]!, cfList:ArrayList[Fragment!]!, combs:Rail[Int], nTerms:Int, sign:Int) : Boolean {
       var fInter:Fragment! = fragList.get(combs(0));       

       for(var i:Int=1; i<=nTerms; i++) {
           fInter = fInter.intersection(fragList.get(combs(i)));

           if (fInter.getNumberOfTrueAtoms() == 0) return false;
       } // end for

       addFragmentToList(fragList, cfList, fInter as Fragment!, sign);
       return true;
   }

   private def addFragmentToList(fragList:ArrayList[Fragment!]!, cfList:ArrayList[Fragment!]!, fragment:Fragment!, sign:Int) {
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

