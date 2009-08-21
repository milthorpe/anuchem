package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a N-dimentional Vector
 * (Initial DRAFT)
 * 
 * @author V.Ganesh
 */
public class Vector { 
    val vec:Array[Double]{rank==1}; 
    val region:Region{rank==1};
    val distribution:Dist{rank==1};

    /**
     * Construct a Vector of dimention N, with default block distribution
     */
    public def this(siz:Int) { 
        region = 0..siz;
        distribution = Dist.makeBlock(region, 0);
        vec = Array.make[Double](distribution);
    }

    /**
     * Construct a Vector from a Matric
     */
    public def this(mat:Matrix) {
        this(mat.getRowCount()*mat.getRowCount());        
         
        var ii:Int = 0;
        val m = mat.getMatrix();
        val N = mat.getRowCount();

        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++)
           for(var j:Int=0; j<N; j++)
              vec(ii++) =  m(i, j);
    }

    /**
     * Construct a Vector of dimention N, with a custom distribution
     */
    public def this(dist:Dist{rank==1}) {
        region = dist.region;
        distribution = dist;
        vec = Array.make[Double](distribution);
    }
 
    /**
     * The size of this matrix
     */
    public def getSize() : Int = vec.region.max(0);
 
    /**
     * Actual data stored
     */
    public def getVector() : Array[Double]{rank==1} = vec;

    /**
     * Get associated region
     */ 
    public def region() : Region{rank==1} = region;
 
    /**
     * Get associated distribution
     */
    public def dist() : Dist{rank==1} = distribution;

    /**
     * make this Vector a null vector
     */
    public def makeZero() : void {
        finish ateach((i) in vec.dist) vec(i) = 0;
    }

    /**
     * the dot product
     */
    public def dot(b:Vector) : Double {
        var res:Double = 0.0;

        // TODO: this is introduced till a fix for XTENLANG-508 is there
        // finish ateach((i) in vec.dist) res += vec(i) * b.vec(i);
        for((i) in vec.dist) res += vec(i) * b.vec(i);

        return res;
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val N = getSize();
        val res = new Vector(N);

        finish ateach((i) in vec.dist) res.vec(i) = vec(i) + b.vec(i);

        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val N = getSize();
        val res = new Vector(N);

        finish ateach((i) in vec.dist) res.vec(i) = vec(i) - b.vec(i);

        return res;
    }

    /**
     * magnitude of this vector
     */
    public def magnitude() : Double {
        var length:Double = 0.0;
        val N = getSize();
        
        finish ateach((i) in vec.dist) length += vec(i) * vec(i);
        
        return Math.sqrt(length);
    }

    /**
     * return a normalized form of this vector
     */
    public def normalize() : Vector {
        val mag = magnitude();
        val N = getSize();
        
        val n = new Vector(N);
       
        finish ateach((i) in vec.dist) n.vec(i) = vec(i) / mag;
        
        return n;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val N = getSize();
        val n = new Vector(N);
        
        finish ateach((i) in vec.dist) n.vec(i) = -vec(i);
        
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {        
        val N = getSize();
        val res = new Vector(N);
        
        finish ateach((i) in vec.dist) res.vec(i) = vec(i) * k;
        
        return res;
    }

    public def toString() : String {
         var str : String = "";
         val N = getSize();

         for(var i:Int=0; i<N; i++) {
            str += vec(i) + " ";
         } // end for

         return str;
    }

}

