/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2013.
 * (C) Copyright Josh Milthorpe 2011.
 */
package x10x.matrix;

import x10.array.Array_2;
import x10.compiler.Inline;

import x10x.vector.Vector;
import x10x.xla.JacobiDiagonalizer;

/**
 * This class represents an (NxM)  Matrix.
 * TODO distributed matrices
 *
 * @author V.Ganesh, milthorpe
 */
public class Matrix { 
    protected val mat:Array_2[Double];

    /**
     * Make instance of Matrix class 
     * 
     * @param siz the size of this matrix
     */
    public def this(siz:Long) {
        mat = new Array_2[Double](siz-1,siz-1);
    }

    /**
     * Make instance of Matrix class
     *
     * @param siz the size of this matrix
     */
    public def this(row:Long, col:Long) {
        mat = new Array_2[Double](row-1, col-1);
    }

    /**
     * Construct a Matrix as a copy of an existing Matrix
     */
    public def this(source : Matrix) {
        mat = new Array_2[Double](source.getMatrix());
    }

    public def getMatrix() = mat;
    public def getRowCount() = mat.numElems_1;
    public def getColCount() = mat.numElems_2;

    public @Inline operator this(i0:Long, i1:Long) : Double {
        return mat(i0,i1);
    }

    public @Inline operator this(i0:Long, i1:Long)=(v:Double) : Double {
        mat(i0,i1) = v;
        return v;
    }

    /**
     * Perform symmetric orthogonalization of this matrix
     */
    public def symmetricOrthogonalization() : Matrix { 
        val diag = new JacobiDiagonalizer();
        diag.diagonalize(this);
        
        val rowCount      = getRowCount(); 
        val eigenValues   = diag.getEigenValues().getVector();
        val eigenVectors  = diag.getEigenVectors();
        val sHalf         = new Matrix(rowCount);

        for (i in 0..(mat.numElems_1-1)) {
            sHalf.mat(i,i) = 1.0 / Math.sqrt(eigenValues(i));
        }

        return (sHalf.similarityTransformT(eigenVectors)); 
    }

    /**
     * Make the current Matrix as Identity
     */
    public def makeIdentity() : void {
        mat.clear();
        for (i in 0..(mat.numElems_1-1)) {
            mat(i,i) = 1.0;
        }
    }

    /**
     * Fill the current Matrix with zero
     */
    public def makeZero() : void {
        mat.clear();
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
    * Multiply two matrices
    * @return C = this . B
    */
    public def mul(b:Matrix) : Matrix {
        val N   = getRowCount();
        val N1  = b.getRowCount();
        val M   = b.getColCount();
        val res = new Matrix(N, M);

        for([i,j] in res.mat.indices()) {
            var cij:Double = 0.0;
            for(var k:Long=0; k<N1; k++) {
               cij += mat(i, k) * b.mat(k, j);
            }
            res.mat(i, j) = cij;
        }

        return res;
    }

   /**
    * Multiply two matrices
    * sets this = A . B
    */
    public def mulInPlace(a:Matrix, b:Matrix) {
        val N1  = b.getRowCount();

        for([i,j] in mat.indices()) {
            var cij:Double = 0.0;
            for(var k:Long=0; k<N1; k++) {
               cij += a.mat(i, k) * b.mat(k, j);
            }
            mat(i, j) = cij;
        }
    }

    /**
     * Scale each element of this matrix by fac
     */
    public def mul(fac:Double) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();

        val res = new Matrix(N, M);
        for([i,j] in mat.indices()) {
            mat(i,j) *= fac;
        }
        //mat.map[Double](mat, (a : Double) => a * fac);
        return res;
    }

    /**
     * Add two matrices
     * @return this + X
     */
    public def add(x:Matrix) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(N, M);
        for([i,j] in mat.indices()) {
            res.mat(i,j) = mat(i,j) + x.mat(i,j);
        }
        //mat.map[Double,Double](res.mat, x.mat, (a : Double, b : Double) => a + b);
        return res;
    }

    /**
     * set this = this + X
     */
    public def addInPlace(x:Matrix) {
        for([i,j] in mat.indices()) {
            mat(i,j) = mat(i,j) + x.mat(i,j);
        }
        //mat.map[Double,Double](mat, x.mat, (a : Double, b : Double) => a + b);
    }

    /**
     * set this = X + Y
     */
    public def addInPlace(x:Matrix, y:Matrix)  {
        for([i,j] in mat.indices()) {
            mat(i,j) = x.mat(i,j) + y.mat(i,j);
        }
        //x.mat.map[Double,Double](mat, y.mat, (a : Double, b : Double) => a + b);
    }

    /**
     * Subtract two matrices: this - X
     */
    public def sub(x:Matrix) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(N, M);
        for([i,j] in mat.indices()) {
            res.mat(i,j) = mat(i,j) - x.mat(i,j);
        }
        //mat.map[Double,Double](res.mat, x.mat, (a : Double, b : Double) => a - b);
        return res;
    }

    /**
     * Find transpose of this matrix
     */
    public def transpose() : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(M, N);

        for([i,j] in res.mat.indices()) {
            res.mat(i, j) = mat(j, i);
        }
 
        return res;
    }

    /**
     * Find trace of this matrix.
     */
    public def trace() : Double {
        var tr : Double = 0.0;

        for (i in 0..(mat.numElems_1-1)) {
            tr += mat(i, i);
        }

        return tr;
    }

    /**
     * absolute sum of off-diagonal elements
     */
    public def sumOffDiagonal() : Double {
        var sum : Double = 0.0;
        val N = getRowCount();

        for (i in 0..(mat.numElems_1-1)) {
            for (j in (i+1)..(mat.numElems_2-1)) {
                sum += Math.abs(mat(i, j));
            }
        }

        return sum; 
    }

    public def getRowVector(rowIdx:Long) : Vector {
       val N = getColCount();

       val vec = new Vector(N);
       val vecVal = vec.getVector();

       for(var i:Long=0; i<N; i++)
           vecVal(i) = mat(rowIdx, i);

       return vec;
    }

    public def isUpperTriangular() : Boolean {
        if (getRowCount() != getColCount()) return false;

        val N = getRowCount();

        for(var i:Long=0; i<N; i++) {
            for(var j:Long=0; j<i; j++) {
                if (mat(i, j) != 0.0) {
                    return false;
                }
            }
        }

        return true;
    }

    public def isSingular(p:Long, row:Rail[Long]) : Boolean {
        val N = getColCount();
        for(i in p..N) {
            if (mat(row(p), i) != 0.0 && Math.abs(mat(row(p),i))>1e-15) {
                return false;
            }
        }

        return true;
    }

    public def toString() : String { 
         var str : String = "";
         val N = getRowCount();
         val M = getColCount();
        
         for(var i:Long=0; i<N; i++) {
            for(var j:Long=0; j<M; j++)  
               str += "" + mat(i, j) + " ";
            str += "\n";
         }

         return str;
    }

   /**
    * Inverse of this matrix
    */
    //@Incomplete public def inverse() : Matrix;

    /**
     * The determinant of this Matrix
     */
    //@Incomplete public def determinant() : Double;

    /**
     * Get a row of this matrix
     */
    //@Incomplete public def getRow(row:Long) : Vector;

    /**
     * Get a column of this matrix
     */
    //@Incomplete public def getColumn(col:Long) : Vector;
}

