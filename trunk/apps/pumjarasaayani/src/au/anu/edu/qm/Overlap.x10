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
    public def this(siz:Int) {
       super(siz);
    }

    var sHalf:Matrix = null;
    public def getSHalf() : Matrix {
       if (sHalf == null) sHalf = this.symmetricOrthogonalization();

       return sHalf;
    }    
}

