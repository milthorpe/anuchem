package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

/**
 * This class represents a toplevel interface for Matrix diagonalizers 
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public interface Diagonalizer {
    public def diagonalize(matrix:Matrix) : void;
    public def getEigenValues()  : Vector;
    public def getEigenVectors() : Matrix;
}

