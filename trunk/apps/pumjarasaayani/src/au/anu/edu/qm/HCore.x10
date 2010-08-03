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

import x10x.matrix.Matrix;

/**
 * HCore.x10
 *
 * The HCore matrix in HF calculation
 *
 * @author: V.Ganesh
 */
public class HCore extends Matrix { 
    public def this(n:Int) {
        super(n);
    }
}

