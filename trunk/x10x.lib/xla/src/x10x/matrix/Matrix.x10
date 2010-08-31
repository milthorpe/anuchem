package x10x.matrix;

import x10x.xla.*;
import x10x.vector.Vector;

import x10.array.Array;

/**
 * This class represents an (NxM)  Matrix.
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public class Matrix { 
    global val mat:Array[Double](2){rect, self.at(this)};

    global val region:Region{rect,rank==2};
    global val distribution:Dist{rect,rank==2}; // TODO actually distribute Matrix!

    /**
     * Make instance of Matrix class 
     * 
     * @param siz the size of this matrix
     */
    public def this(siz:Int) {
        region       = [0..(siz-1), 0..(siz-1)];
        // distribution = Dist.makeBlock(region, 1);
        distribution = Dist.makeConstant(region);
        mat          = new Array[Double](region);
    }

    /**
     * Make instance of Matrix class
     *
     * @param siz the size of this matrix
     */
    public def this(row:Int, col:Int) {
        region       = [0..(row-1), 0..(col-1)];
        // distribution = Dist.makeBlock(region, 1);
        distribution = Dist.makeConstant(region);
        mat          = new Array[Double](region);
    }

    /**
     * Construct a Matrix with a custom distribution
     */
    public def this(dist:Dist{rank==2}) {
        distribution = dist;
        region       = distribution.region;
        mat          = new Array[Double](region);
    }

    /**
     * Construct a Matrix as a copy of an existing Matrix (may be remote)
     */
    public def this(source : Matrix) {
        val sourceMat = source.getMatrix();
        distribution = source.distribution;
        region       = sourceMat.region;
        mat          = new Array[Double](sourceMat.region);
        finish {mat.copyFrom(sourceMat);}
    }

    public global def region() = region;
    public global def dist() = distribution;

    public global def getMatrix() = mat;
    public global def getRowCount() = mat.region.max(0)+1;
    public global def getColCount() = mat.region.max(1)+1;

    /**
     * Perform symmetric orthogonalization of this matrix
     */
    public def symmetricOrthogonalization() : Matrix! { 
        val diag = new JacobiDiagonalizer();
       
        diag.diagonalize(this);
        
        val rowCount      = getRowCount(); 
        val eigenValues   = diag.getEigenValues().getVector();
        val eigenVectors  = diag.getEigenVectors();
        val sHalf         = new Matrix(rowCount);
        
        sHalf.makeIdentity();

        val sqrtEVal = Rail.make[Double](rowCount);

        finish foreach((i,j) in sHalf.mat.region) {
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
        finish foreach((i,j) in mat)
           if (i == j) mat(i, j) = 1.0;
           else        mat(i, j) = 0.0;
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
    public def similarityTransformT(x:Matrix{self.at(this)}) : Matrix{self.at(this)} {
        return x.transpose().mul(this).mul(x);
    }

    /**
     * Perform a similarity transform: X . this . X'
     */
    public def similarityTransform(x:Matrix{self.at(this)}) : Matrix{self.at(this)} {
        return x.mul(this).mul(x.transpose());
    }

    /**
     * Multiply two matrices: this . X
     */
    public def mul(x:Matrix!) : Matrix! {
         val N   = getRowCount();
         val N1  = x.getRowCount();
         val M   = x.getColCount();
         val res = new Matrix(N, M);

         for((i, j) in res.mat) {
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
    public def mul(fac:Double) : Matrix! {
         val N   = getRowCount();
         val M   = getColCount();

         val res = new Matrix(N, M);
      
         finish foreach((i, j) in res.mat)
            res.mat(i, j) = mat(i, j) * fac;

         return res;
    }

    /**
     * Add two matrices: this + X
     */
    public def add(x:Matrix!) : Matrix! {
         val N   = getRowCount();
         val M   = getColCount();
         val res = new Matrix(N, M);

         finish foreach((i, j) in res.mat)
            res.mat(i, j) = mat(i, j) + x.mat(i, j);

         return res;
    }

    /**
     * Subtract two matrices: this - X
     */
    public def sub(x:Matrix!) : Matrix! {
         val N   = getRowCount();
         val M   = getColCount();
         val res = new Matrix(N, M);

         finish foreach((i, j) in res.mat)
            res.mat(i, j) = mat(i, j) - x.mat(i, j);

         return res;
    }

    /**
     * Find transpose of this matrix
     */
    public def transpose() : Matrix! {
         val N   = getRowCount();
         val M   = getColCount();
         val res = new Matrix(M, N);

         finish foreach((i,j) in res.mat) {
             res.mat(i, j) = mat(j, i);
         }
 
         return res;
    }

    /**
     * Find trace of this matrix.
     */
    public def trace() : Double {
         val N = getRowCount();
         val tr = Rail.make[Double](1);
         
         tr(0) = 0.0;
         finish foreach((i,j) in mat) {
                 if (i==j) {
                     atomic tr(0) += mat(i, i);
                 }
         }
        
         return tr(0);
    }

    /**
     * absolute sum of off-diagonal elements
     */
    public def sumOffDiagonal() : Double {
       // sum off diagonal elements of A
       val sum = Rail.make[Double](1);
       val N = getRowCount();

       sum(0) = 0.0;
       finish foreach((i,j) in mat) {
                if (i!=j && j>i) {
                    atomic sum(0) += Math.abs(mat(i, j));
                }
       }

       return sum(0); 
    }

    public def getRowVector(rowIdx:Int) : Vector! {
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

    public def isSingular(p:Int, row:Rail[Int]!) : Boolean {
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

    public def getValRail() {
        val N = getRowCount();
        val M = getColCount();
        val r = Rail.make[Double](N*M);
        var i:Int, j:Int, ii:Int;

        ii = 0;
        for(i=0; i<N; i++)
           for(j=0; j<M; j++)
              r(ii++) = mat(i,j);

        return ValRail.make(r);
    }

    public global safe def toString() : String { 
         var str : String = "";
         val N = getRowCount();
         val M = getColCount();
        
         // TODO: 
         for(var i:Int=0; i<N; i++) {
            for(var j:Int=0; j<M; j++)  
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

