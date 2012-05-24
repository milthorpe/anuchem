/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.matrix.DenseMatrix;

/**
 * Represents the overlap matrix
 *
 * @author: V.Ganesh
 */
public class Overlap extends DenseMatrix{self.M==self.N} {
    var sHalf:DenseMatrix{self.M==this.M,self.N==this.N} = null;

    public def this(n:Int):Overlap{self.M==n,self.N==n} {
        super(n, n);
    }

    public def getSHalf():DenseMatrix(this.N, this.N) {
       if (sHalf == null) sHalf = GMLDiagonalizer.symmetricOrthogonalization(this);

       return sHalf;
    }    
}

