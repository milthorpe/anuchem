package x10x.matrix;

import x10x.xla.*;
import x10x.vector.Vector;

/**
 * This class represents an (NxM)  Matrix.
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public class Matrix { 
    var mat : Array[Double]{rank==2};

    /**
     * Construct an NxN matrix
     *
     * @param siz the size of this matrix
     */
    public def this(siz:Int) {
        mat = Array.make[Double]([0..siz, 0..siz]);
    }

    /**
     * Construct an NxM matrix
     *
     * @param row number of rows in this matrix
     * @param col number of rows in this matrix
     */
    public def this(row:Int, col:Int) {
        mat = Array.make[Double]([0..row, 0..col]);
    }

    public def getMatrix() : Array[Double]{rank==2} = mat;
    public def getRowCount() : Int = mat.region.max(0);
    public def getColCount() : Int = mat.region.max(1);

    /**
     * Perform symmetric orthogonalization of this matrix
     */
    public def symmetricOrthogonalization() : Matrix { 
        val diag = new JacobiDiagonalizer();
        
        diag.diagonalize(this);
        
        val rowCount     = getRowCount(); 
        val eigenValues  = diag.getEigenValues();
        val eigenVectors = diag.getEigenVectors();
        val sHalf        = new Matrix(rowCount);
        
        sHalf.makeIdentity();

        // TODO : x10 parallel
        for(var i:Int=0; i<rowCount; i++) 
            sHalf.mat(i, i) /= Math.sqrt(eigenValues(i));
        
        return (sHalf.similarityTransformT(eigenVectors)); 
    }
  
    /**
     * Make the current Matrix as Identity
     */
    public def makeIdentity() : void {
        val N = getRowCount();

        // TODO : x10 parallel
        for(var i:Int=0; i<N; i++) 
           for(var j:Int=0; j<N; j++)
              if (i == j) mat(i, j) = 1.0;
              else        mat(i, j) = 0.0;
    }

    /**
     * Perform a similarity transform: X' . this . X
     */
    public def similarityTransformT(x:Matrix) : Matrix {
        return x.transpose().mul(this).mul(x);
    }

    /**
     * Perform a similarity transform: X . this . X'
     */
    public def similarityTransform(x:Matrix) : Matrix {
        return x.mul(this).mul(x.transpose());
    }

    /**
     * Multiply two matrices: this . X
     */
    public def mul(x:Matrix) : Matrix {
         val N  = getRowCount();
         val N1 = x.getRowCount();
         val M  = x.getColCount();
         val res = new Matrix(N, M);
         var cij:Double;

         // TODO : x10 parallel
         for(var i:Int=0; i<N; i++) 
            for(var j:Int=0; j<M; j++) {
               cij = 0.0;
               for(var k:Int=0; k<N1; k++) 
                  cij += mat(i, k) * x.mat(k, j);
               res.mat(i, j) = cij;
            } // end for  

         return res;
    }

    /**
     * Add two matrices: this + X
     */
    public def add(x:Matrix) : Matrix {
         val N = getRowCount();
         val M = getColCount();
         val res = new Matrix(N, M);

         // TODO : x10 parallel
         for(var i:Int=0; i<N; i++) 
            for(var j:Int=0; j<M; j++) 
               res.mat(i, j) = mat(i, j) + x.mat(i, j);

         return res;
    }

    /**
     * Subtract two matrices: this - X
     */
    public def sub(x:Matrix) : Matrix {
         val N = getRowCount();
         val M = getColCount();
         val res = new Matrix(N, M);

         // TODO : x10 parallel
         for(var i:Int=0; i<N; i++)
            for(var j:Int=0; j<M; j++)
               res.mat(i, j) = mat(i, j) - x.mat(i, j);

         return res;
    }

    /**
     * Find transpose of this matrix
     */
    public def transpose() : Matrix {
         val N = getRowCount();
         val M = getColCount();
         val res = new Matrix(M,N);

         // TODO : x10 parallel
         for(var i:Int=0; i<N; i++) {
            for(var j:Int=0; j<M; j++) {
               res.mat(j, i) = mat(i, j);
            } // end for
         } // end for i 
 
         return res;
    }

    /**
     * Find trace of this matrix.
     */
    public def trace() : Double {
         val N = getRowCount();
         var tr:Double = 0.0;
        
         // TODO : x10 parallel
         for(var i:Int=0; i<N; i++) 
            tr += mat(i, i);
        
         return tr;
    }

    public def toString() : String { 
         var str : String = "";
         val N = getRowCount();
         
         for(var i:Int=0; i<N; i++) {
            for(var j:Int=0; j<N; j++)  
               str += mat(i, j) + " ";
            str += "\n";
         } // end for

         return str;
    }

   /**
    * Inverse of this matrix
    */
    public incomplete def inverse() : Matrix;

    /**
     * The determinant of this Matrix
     */
    public incomplete def determinant() : Double;

    /**
     * Get a row of this matrix
     */
    public incomplete def getRow(row:Int) : Vector;

    /**
     * Get a column of this matrix
     */
    public incomplete def getColumn(col:Int) : Vector;
}

