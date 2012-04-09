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

import x10.matrix.DenseMatrix;

/**
 * The Fock matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Fock extends DenseMatrix {
    public def this(n:Int):Fock{self.M==n,self.N==n} {
        super(n,n);
    }
  
    public def compute(hCore:HCore{self.M==this.M,self.N==this.N}, gMatrix:DenseMatrix{self.M==this.M,self.N==this.N}):void {
        hCore.copyTo(this);
        super.cellAdd(gMatrix);
    }
}

