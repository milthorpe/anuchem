/**
 * Utility.x10
 *
 * A few utility functions
 *
 */

package au.anu.edu.qm;

import x10.util.*;

public class Utility {

    public def this() { }

    /**
     * split a string and return the result as an ArrayList
     */
    public static def split(str:String, splitChar:Char) : ArrayList[String] {
        val strs = new ArrayList[String]();

        var ns:String = str;
        while(true) {
           val i = ns.indexOf(splitChar);
           if (i < 0) {
             if(ns.length() != 0) strs.add(ns);
             break;
           } else if (i == 0) {
             ns = ns.substring(i+1, ns.length());
           } else {
             strs.add(ns.substring(0, i));
             ns = ns.substring(i+1, ns.length());
           } // end if
        }

        return strs;
    }
}

