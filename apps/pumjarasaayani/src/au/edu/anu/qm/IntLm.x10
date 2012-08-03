/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.compiler.NonEscaping;
import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * Represents a list of ShellPair
 *
 * @author: T. Limpanuparb, J. Milthorpe
 */

public struct IntLm {
    public val IntLm : Rail[Double];
    public val L : Int;

    public def this(IntLm:Rail[Double], L:Int) {
        this.IntLm=IntLm;
        this.L=L;
    }

}

