package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a N-dimentional Vector
 * (Initial DRAFT)
 * TODO distribute Vector
 * 
 * @author V.Ganesh
 */
public class Vector { 
    val vec:Array[Double](1);
    val region:Region(1){rect};

    /**
     * Construct a Vector of dimention N
     */
    public def this(siz:Int) { 
        region       = [0..(siz-1)];
        // distribution = Dist.makeBlock(region, 0);
        vec          = new Array[Double](region);
    }

    /**
     * Construct a Vector from a Matrix
     */
    public def this(mat:Matrix) {
        this(mat.getRowCount()*mat.getColCount());
         
        var ii:Int = 0;
        val m = mat.getMatrix();

        for([i,j] in m.region)
            vec(ii++) =  m(i, j);
    }
 
    /**
     * The size of this matrix
     */
    public def getSize() = vec.region.max(0)+1;
 
    /**
     * Actual data stored
     */
    public def getVector() = vec;

    /**
     * Get associated region
     */ 
    public def region() = region;

    /**
     * make this Vector a null vector
     */
    public def makeZero() : void {
        vec.fill(0.0);
    }

    /**
     * the dot product
     */
    public def dot(b:Vector) : Double {
        val res = Rail.make[Double](1);

        // TODO: this is introduced till a fix for XTENLANG-508 is there
        // finish foreach([i] in vec.region) res += vec(i) * b.vec(i);
        finish foreach([i] in vec) res(0) += vec(i) * b.vec(i);

        return res(0);
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val N   = getSize();
        val res = new Vector(N);

        finish foreach([i] in vec) res.vec(i) = vec(i) + b.vec(i);

        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val N   = getSize();
        val res = new Vector(N);

        finish foreach([i] in vec) res.vec(i) = vec(i) - b.vec(i);

        return res;
    }

    /**
     * magnitude of this vector
     */
    public def magnitude() : Double {
        var magSquared: Double = 0.0;
        val N = getSize();

        for (var i:Int=0; i<N; i++) 
            magSquared += vec(i)*vec(i);
       
        return Math.sqrt(magSquared);
    }

    /**
     * return a normalized form of this vector
     */
    public def normalize() : Vector {
        val mag = magnitude();
        val N   = getSize();
        
        val n = new Vector(N);
       
        finish foreach([i] in vec) n.vec(i) = vec(i) / mag;
        
        return n;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val N = getSize();
        val n = new Vector(N);
        
        finish foreach([i] in vec) n.vec(i) = -vec(i);
        
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {        
        val N   = getSize();
        val res = new Vector(N);
        
        finish foreach([i] in vec) res.vec(i) = vec(i) * k;
        
        return res;
    }

    /**
     * max norm of this vector
     */
    public def maxNorm() : Double {
        val N   = getSize();
        val res = Rail.make[Double](1);

        finish foreach([i] in vec) res(0) = Math.max(Math.abs(vec(i)), res(0));
        
        return res(0);       
    }

    public safe def toString() : String {
         var str : String = "";
         val N = getSize();

         for(var i:Int=0; i<N; i++) {
            str += vec(i) + " ";
         } // end for

         return str;
    }

}

