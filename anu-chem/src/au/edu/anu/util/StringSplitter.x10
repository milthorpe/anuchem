/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */
package au.edu.anu.util;

import x10.util.ArrayList;

public class StringSplitter {
    public static def splitOnWhitespace(s:String):Rail[String] {
        val words = new ArrayList[String]();
        val end = s.length();
        var offset:Int = 0n;
        var startOffset:Int=0n;
        var word:Boolean = !(s.charAt(0n).isWhitespace());
        while(offset < end) {
            val whitespace = s.charAt(offset).isWhitespace();
            if (word && whitespace) {
                word = false;
                words.add(s.substring(startOffset, offset));
            } else if (!word && !whitespace) {
                word = true;
                startOffset = offset;
            }
            offset++;
        }
        if (word) words.add(s.substring(startOffset, offset));
        return words.toRail();
    }
}

