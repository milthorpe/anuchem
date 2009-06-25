package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a N-dimentional Vector
 * (Initial DRAFT)
 * 
 * @author V.Ganesh
 */
public class Vector { 
    var vec:Array[Double]{rank==1}; 

    /**
     * Construct a Vector of dimention N
     */
    public def this(siz:Int) { 
        vec = Array.make[Double]([0..siz]);
    } 

    /**
     * Construct a Vector from a Matric
     */
    public def this(mat:Matrix) { 
        val siz = mat.getRowCount();
        vec = Array.make[Double]([0..siz*siz]);
        
        var ii:Int = 0;
        val m = mat.getMatrix();

        // TODO : x10 - parallel
        for(var i:Int=0; i<siz; i++)
           for(var j:Int=0; j<siz; j++)
              vec(ii++) =  m(i, j);
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
     * make this Vector a null vector
     */
    public def makeZero() : void {
        val N = getSize();

        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) vec(i) = 0;
    }

    /**
     * the dot product
     */
    public def dot(b:Vector) : Double {
        val N = getSize();
        var res:Double = 0.0;
      
        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) res += vec(i) * b.vec(i);
 
        return res;
    }

    /**
     * add two vectors: this + X 
     */
    public def add(x:Vector) : Vector {
        val N = getSize();
        val res = new Vector(N);

        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) res.vec(i) = vec(i) + b.vec(i);

        return res;
    }

    /**
     * subtract two vectors: this - X
     */
    public def sub(x:Vector) : Vector {
        val N = getSize();
        val res = new Vector(N);

        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) res.vec(i) = vec(i) - b.vec(i);

        return res;
    }

    /**
     * magnitude of this vector
     */
    public def magnitude() : Double {
        var length:Double = 0.0;
        val N = getSize();
        
        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) length += vec(i) * vec(i);
        
        return Math.sqrt(length);
    }

    /**
     * return a normalized form of this vector
     */
    public def normalize() : Vector {
        val mag = magnitude();
        val N = getSize();
        
        val n = new Vector(N);
       
        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) n.vec(i) = vec(i) / mag;
        
        return n;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {      
        val n = new Vector(vector.length);
        val N = getSize();
        
        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++) n.vec(i) = -vec(i);
        
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {        
        val N = getSize();
        val result = new Vector(N);
        
        for(var i:Int=0; i<N; i++) result.vec(i) = this.vec(i) * k;
        
        return result;
    }
}

