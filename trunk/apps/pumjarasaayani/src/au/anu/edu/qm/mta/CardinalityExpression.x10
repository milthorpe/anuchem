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
   public def this() {
   }

   /** add cardinality fragments */
   public def addCardinalityFragments(fragList:ArrayList[Fragment]!) {
       val noOfFragments = fragList.size();
       val combs = Rail.make[Int](noOfFragments, (Int)=>0);

       for(var i:Int=0; i<noOfFragments; i++) {
          for(var j:Int=i+1; j<noOfFragments; j++) {
             combs(0) = i; combs(1) = j;

             val sign = getSignOfTerm(2); 

              
          } // end for
       } // end for
   }

   private def getSignOfTerm(n:Int) : Int {
       if (n%2 == 0) return -1;
       else          return 1;
   }

   private def computeIntersections(combs:Rail[Int], nTerms:Int, sign:Int) : Boolean {
       if (nTerms == 0) return false;
 
       var fInter:Fragment = combs(0);       

       for(i:Int=1; i<nTerms; i++) {
           fInter = fInter.intersection(combs(i));
       } // end for

       return true;
   }
}

