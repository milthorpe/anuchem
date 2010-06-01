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

package au.anu.edu.qm.mta;

import x10.util.ArrayList;
import x10.util.ValRailBuilder;

public class CardinalityExpression {

   private var mergedIntoPreviousTerm:Boolean;

   public def this() {
   }

   /** add cardinality fragments */
   public def addCardinalityFragments(fragList:ArrayList[Fragment]!) {
       val noOfFragments = fragList.size();
       val combs = Rail.make[Int](noOfFragments, (Int)=>0);
       val cfList = new ArrayList[Fragment]();

       for(var i:Int=0; i<noOfFragments; i++) {
          for(var j:Int=i+1; j<noOfFragments; j++) {
             combs(0) = i; combs(1) = j;

             val sign = getSignOfTerm(2); 

             if (!computeIntersections(fragList, cfList, combs, 2, sign)) continue;

             // TODO:              
          } // end for
       } // end for

       for(cf in cfList) fragList.add(cf);
   }

   private def getSignOfTerm(n:Int) : Int {
       if (n%2 == 0) return -1;
       else          return 1;
   }

   private def computeIntersections(fragList:ArrayList[Fragment]!, cfList:ArrayList[Fragment]!, combs:Rail[Int], nTerms:Int, sign:Int) : Boolean {
       if (nTerms == 0) return false;
 
       var fInter:Fragment = fragList.get(combs(0));       

       for(var i:Int=1; i<nTerms; i++) {
           fInter = fInter.intersection(fragList.get(combs(i)) as Fragment!);
       } // end for

       if (fInter.getNumberOfAtoms() == 0) return false;
 
       addFragmentToList(cfList, fInter as Fragment!, sign);
       return true;
   }

   private def addFragmentToList(cfList:ArrayList[Fragment]!, fragment:Fragment!, sign:Int) {
       mergedIntoPreviousTerm = false;

       for(frag in cfList) {
          val nCommonAtoms = frag.intersection(fragment).getNumberOfAtoms();

          if (nCommonAtoms == frag.getNumberOfTrueAtoms()) {
             mergedIntoPreviousTerm = true;
             frag.cardinalitySign(frag.cardinalitySign() + sign);
          } // end if
       } // end for
   }
}

