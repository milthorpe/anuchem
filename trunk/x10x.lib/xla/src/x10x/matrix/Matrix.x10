package x10x.matrix;

import x10x.xla.*;
import x10x.vector.Vector;

import x10.array.Array;
import x10.compiler.Header;
import x10.compiler.Inline;

/**
 * This class represents an (NxM)  Matrix.
 * (Initial DRAFT)
 * TODO distribute Matrix
 *
 * @author V.Ganesh
 */
public class Matrix { 
    protected val mat:Array[Double]{rect,rank==2};
    val region:Region{rect,rank==2};

    /**
     * Make instance of Matrix class 
     * 
     * @param siz the size of this matrix
     */
    public def this(siz:Int) {
        region       = (0..(siz-1)) * (0..(siz-1));
        mat          = new Array[Double](region);
    }

    /**
     * Make instance of Matrix class
     *
     * @param siz the size of this matrix
     */
    public def this(row:Int, col:Int) {
        region       = (0..(row-1)) * (0..(col-1));
        mat          = new Array[Double](region);
    }

    /**
     * Construct a Matrix as a copy of an existing Matrix
     */
    public def this(source : Matrix) {
        val sourceMat = source.getMatrix();
        mat = sourceMat;
        region       = sourceMat.region;
        //mat          = new Array[Double](sourceMat.region);
        //finish {mat.asyncCopy(sourceMat);}
    }

    public def region() = region;

    public def getMatrix() = mat;
    public def getRowCount() = mat.region.max(0)+1;
    public def getColCount() = mat.region.max(1)+1;

    public @Inline operator this(i0:int, i1:int) : Double {
        return mat(i0,i1);
    }

    public @Inline operator this(i0:int, i1:int)=(v:Double) : Double {
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
        
        sHalf.makeIdentity();

        finish for([i,j] in sHalf.mat.region) async {
             if (i==j) {
                sHalf.mat(i,i) /= Math.sqrt(eigenValues(i));
             }
        }

        return (sHalf.similarityTransformT(eigenVectors)); 
    }

    /**
     * Make the current Matrix as Identity
     */
    public def makeIdentity() : void {
        mat.fill(0.0);
        for ([i] in 0..mat.region.max(0)) {
            mat(i,i) = 1.0;
        }
    }

    /**
     * Fill the current Matrix with zero
     */
    public def makeZero() : void {
        mat.fill(0.0);
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
         val N   = getRowCount();
         val N1  = x.getRowCount();
         val M   = x.getColCount();
         val res = new Matrix(N, M);

         for([i,j] in res.mat) {
            var cij:Double = 0.0;
            for(var k:Int=0; k<N1; k++) {
               cij += mat(i, k) * x.mat(k, j);
            }
            res.mat(i, j) = cij;
         } // end for  

         return res;
    }

    /**
     * Scale each element of this matrix by fac
     */
    public def mul(fac:Double) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();

        val res = new Matrix(N, M);

        val mat = this.mat; // TODO this should not be required XTENLANG-1913
        finish for([i,j] in res.mat) async
            res.mat(i, j) = mat(i, j) * fac;

        return res;
    }

    /**
     * Add two matrices: this + X
     */
    public def add(x:Matrix) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(N, M);

        val mat = this.mat; // TODO this should not be required XTENLANG-1913
        finish for([i,j] in res.mat) async
            res.mat(i, j) = mat(i, j) + x.mat(i, j);

        return res;
    }

    /**
     * this = this + X
     */
    public def addInPlace(x:Matrix)  {
        val N   = getRowCount();
        val M   = getColCount();

        val mat = this.mat; // TODO this should not be required XTENLANG-1913
        finish for([i,j] in mat) async
            mat(i, j) = mat(i, j) + x.mat(i, j);
    }

    /**
     * Subtract two matrices: this - X
     */
    public def sub(x:Matrix) : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(N, M);

        val mat = this.mat; // TODO this should not be required XTENLANG-1913
        finish for([i,j] in res.mat) async
            res.mat(i, j) = mat(i, j) - x.mat(i, j);

        return res;
    }

    /**
     * Find transpose of this matrix
     */
    public def transpose() : Matrix {
        val N   = getRowCount();
        val M   = getColCount();
        val res = new Matrix(M, N);

        val mat = this.mat; // TODO this should not be required XTENLANG-1913
        finish for([i,j] in res.mat) async {
            res.mat(i, j) = mat(j, i);
        }
 
        return res;
    }

    /**
     * Find trace of this matrix.
     */
    public def trace() : Double {
        var tr : Double = 0.0;

        for ([i] in 0..mat.region.max(0)) {
            tr += mat(i, i);
        }

        return tr;
    }

    /**
     * absolute sum of off-diagonal elements
     */
    public def sumOffDiagonal() : Double {
       // sum off diagonal elements of A
       val sum = new Array[Double](1);
       val N = getRowCount();

       sum(0) = 0.0;
        val mat = this.mat; // TODO this should not be required XTENLANG-1913
       finish for([i,j] in mat) async {
                if (i!=j && j>i) {
                    atomic sum(0) += Math.abs(mat(i, j));
                }
       }

       return sum(0); 
    }

    public def getRowVector(rowIdx:Int) : Vector {
       val N = getColCount();

       val vec = new Vector(N);
       val vecVal = vec.getVector();

       // TODO:
       for(var i:Int=0; i<N; i++)
           vecVal(i) = mat(rowIdx, i);

       return vec;
    }

    public def isUpperTriangular() : Boolean {
        if (getRowCount() != getColCount()) return false;

        val N = getRowCount();
        var i:Int, j:Int;

        for(i=0; i<N; i++) {
            for(j=0; j<i; j++) {
                if (mat(i, j) != 0) {
                    return false;
                } // end if
            } // end for
        } // end for

        return true;
    }

    public def isSingular(p:Int, row:Array[Int](1){rail}) : Boolean {
        val N = getColCount();
        var i:Int;

        Console.OUT.println("isSingular: " + p + " " + row(p) + ", " + N);
        Console.OUT.println("isSingular: " + this);

        for(i=p; i<=N; i++) {
            if (mat(row(p), i) != 0.0 && Math.abs(mat(row(p),i))>1e-15) {
                return false;
            } // end if
        } // end for

        return true;
    }

    public def toString() : String { 
         var str : String = "";
         val N = getRowCount();
         val M = getColCount();
        
         // TODO: 
         for(var i:Int=0; i<N; i++) {
            for(var j:Int=0; j<M; j++)  
               str += "" + mat(i, j) + " ";
            str += "\n";
         } // end for


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
    //@Incomplete public def getRow(row:Int) : Vector;

    /**
     * Get a column of this matrix
     */
    //@Incomplete public def getColumn(col:Int) : Vector;
}

