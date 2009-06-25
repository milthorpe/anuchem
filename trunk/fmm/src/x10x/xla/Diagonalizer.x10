package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a toplevel interface for Matrix diagonalizers 
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public interface Diagonalizer {
    public def diagonalize(matrix:Matrix) : void;
    public def getEigenValues()  : Array[Double]{rank==1};
    public def getEigenVectors() : Matrix;
}

