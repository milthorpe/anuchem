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
    global val mat:Array[Double]{rank==2, self.at(this)};
    global val region:Region{rank==2, self.at(this)};
    global val distribution:Dist{rank==2, self.at(this)};

    /**
     * Make instance of Matrix class 
     * 
     * @param siz the size of this matrix
     */
    public def this(siz:Int) {
        region       = [0..(siz-1), 0..(siz-1)];
        // distribution = Dist.makeBlock(region, 1) as Dist{rank==2, self.at(this)};
        distribution = Dist.makeConstant(region) as Dist{rank==2, self.at(this)};
        mat          = Array.make[Double](distribution) as Array[Double]{rank==2, self.at(this)};
    }

    /**
     * Make instance of Matrix class
     *
     * @param siz the size of this matrix
     */
    public def this(row:Int, col:Int) {
        region       = [0..(row-1), 0..(col-1)];
        // distribution = Dist.makeBlock(region, 1) as Dist{rank==2, self.at(this)};
        distribution = Dist.makeConstant(region) as Dist{rank==2, self.at(this)};
        mat          = Array.make[Double](distribution) as Array[Double]{rank==2, self.at(this)};
    }

    /**
     * Construct a Matrix with a custom distribution
     */
    public def this(dist:Dist{rank==2}) {
        distribution = dist as Dist{rank==2, self.at(this)};
        region       = distribution.region as Region{rank==2, self.at(this)};
        mat          = Array.make[Double](distribution) as Array[Double]{rank==2, self.at(this)};
    }

    /**
     * Get associated region
     */
    public def region() = region;

    /**
     * Get associated distribution
     */
    public def dist() = distribution;

    public def getMatrix() = mat;
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
        val sHalf         = new Matrix(rowCount) as Matrix!;
        
        sHalf.makeIdentity();

        val sqrtEVal = Rail.make[Double](rowCount);

        finish foreach(plc in eigenValues.dist.places()) {
             for(val(i) in eigenValues.dist.get(plc)) {
                 sqrtEVal(i) = at(plc) { return Math.sqrt(eigenValues(i)); };
             }
        }

        finish foreach(plc in sHalf.mat.dist.places()) {
             for(val(i, j) in sHalf.mat.dist.get(plc)) {
                 if (i == j) { 
                    val sqVal = sqrtEVal(i);
                    at(plc) { sHalf.mat(i, i) /= sqVal; };
                 }
             }
        }

        return (sHalf.similarityTransformT(eigenVectors) as Matrix!); 
    }

    /**
     * Make the current Matrix as Identity
     */
    public def makeIdentity() : void {
        finish ateach(val(i,j) in mat.dist)
           if (i == j) mat(i, j) = 1.0;
           else        mat(i, j) = 0.0;
    }

    /**
     * Fill the current Matrix with zero
     */
    public def makeZero() : void {
        finish ateach(val(i,j) in mat.dist)
           mat(i, j) = 0.0;
    }

    /**
     * Perform a similarity transform: X' . this . X
     */
    public def similarityTransformT(x:Matrix{self.at(this)}) : Matrix{self.at(this)} {
        return x.transpose().mul(this).mul(x) as Matrix{self.at(this)};
    }

    /**
     * Perform a similarity transform: X . this . X'
     */
    public def similarityTransform(x:Matrix{self.at(this)}) : Matrix{self.at(this)} {
        return x.mul(this as Matrix{self.at(this)}).mul(x.transpose());
    }

    /**
     * Multiply two matrices: this . X
     */
    public def mul(x:Matrix!) : Matrix! {
         val N   = getRowCount();
         val N1  = x.getRowCount();
         val M   = x.getColCount();
         val res = new Matrix(N, M);

         for(val(i, j) in res.mat.region) {
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
      
         finish foreach(val(i, j) in res.mat.region)
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

         finish foreach(val(i, j) in res.mat.region)
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

         finish foreach(val(i, j) in res.mat.region)
            res.mat(i, j) = mat(i, j) - x.mat(i, j);

         return res;
    }

    /**
     * Find transpose of this matrix
     */
    public def transpose() : Matrix! {
         val N   = getRowCount();
         val M   = getColCount();
         val res = new Matrix(Dist.make([0..(M-1), 0..(N-1)])) as Matrix!;

         var i:Int, j:Int;

         finish foreach(plc in mat.dist.places()) {
             for(val(i, j) in mat.dist.get(plc)) {
                 res.mat(j, i) = at(plc) { return mat(i, j); }; 
             }
         }

         val distRes = new Matrix(M, N) as Matrix!;

         finish foreach(plc in distRes.mat.dist.places()) {
             for(val(i, j) in distRes.mat.dist.get(plc)) {
                val r = res.mat(i, j);
                at(plc) { distRes.mat(i, j) = r; }
             }
         }
 
         return distRes;
    }

    /**
     * Find trace of this matrix.
     */
    public def trace() : Double {
         val N = getRowCount();
         val tr = Rail.make[Double](1);
         
         tr(0) = 0.0;
         finish foreach(plc in mat.dist.places()) {
             for(val(i, j) in mat.dist.get(plc)) {
                 if (i==j) {
                     val trLoc = at(plc) { return mat(i, i); };

                     atomic tr(0) += trLoc;
                 }
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
       finish foreach(plc in mat.dist.places()) {
             for(val(i, j) in mat.dist.get(plc)) {
                 if (i!=j && j>i) {
                     val sumLoc = at(plc) { return Math.abs(mat(i, j)); };
                     
                     atomic sum(0) += sumLoc;
                 }
             }
       }

       return sum(0); 
    }

    public def getRowVector(rowIdx:Int) : Vector! {
       val N = getColCount();

       val vec = new Vector(N) as Vector!;
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

