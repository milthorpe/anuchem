/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2013.
 */
package au.edu.anu.qm;

/**
 * Represents a list of ShellPair
 *
 * @author: T. Limpanuparb, J. Milthorpe
 */

public struct Ylm {
    public val y : Rail[Double];
    public val maxL : Int;

    public def this(y:Rail[Double], maxL:Int) {
        this.y=y;
        this.maxL=maxL;
    }

}

