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

