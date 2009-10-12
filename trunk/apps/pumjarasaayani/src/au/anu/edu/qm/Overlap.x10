/**
 * Overlap.x10
 *
 * Represents the overlap matrix
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;

public class Overlap extends Matrix {
    var sHalf:Matrix{self.at(this)} = null;

    public def getSHalf() : Matrix{self.at(this)} {
       if (sHalf == null) sHalf = this.symmetricOrthogonalization();

       return sHalf;
    }    
}

