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
package au.anu.edu.qm;

import x10.util.*;

/**
 * Utility.x10
 *
 * A few utility functions
 *
 */
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

