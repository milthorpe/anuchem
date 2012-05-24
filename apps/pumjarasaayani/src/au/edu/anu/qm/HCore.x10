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
package au.edu.anu.qm;

import x10.matrix.DenseMatrix;

/**
 * The HCore matrix in HF calculation
 *
 * @author: V.Ganesh
 */
public class HCore extends DenseMatrix{self.M==self.N} { 
    public def this(n:Int):HCore{self.M==n,self.N==n} {
        super(n, n);
    }
}

