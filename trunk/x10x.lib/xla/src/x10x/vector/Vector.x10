package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a N-dimentional Vector
 * (Initial DRAFT)
 * 
 * @author V.Ganesh
 */
public class Vector { 
    global var vec:Array[Double]{rank==1}; 
    global var region:Region{rank==1};
    global var distribution:Dist{rank==1};

    /**
     * Empty constructor
     */
    public def this() { }

    /**
     * Construct a Vector of dimention N, with default block distribution
     */
    public static def make(siz:Int) : Vector { 
        val newVector = new Vector();

        newVector.region       = 0..(siz-1);
        newVector.distribution = Dist.makeBlock(newVector.region, 0);
        newVector.vec          = Array.make[Double](newVector.distribution);

        return newVector;
    }

    /**
     * Construct a Vector of dimention N, with default block distribution
     */
    public static def make(newVector:Vector, siz:Int) : Vector {
        newVector.region       = 0..(siz-1);
        newVector.distribution = Dist.makeBlock(newVector.region, 0);
        newVector.vec          = Array.make[Double](newVector.distribution);

        return newVector;
    }

    /**
     * Construct a Vector from a Matrix
     */
    public def make(mat:Matrix{self.at(this)}) : Vector {
        val N = mat.getRowCount();
        val newVector = Vector.make(N*N);
         
        var ii:Int = 0;
        val m = mat.getMatrix();

        // TODO : x10 - parallel
        for(var i:Int=0; i<N; i++)
           for(var j:Int=0; j<N; j++)
              newVector.vec(ii++) =  m(i, j);

        return newVector;
    }

    /**
     * Construct a Vector of dimention N, with a custom distribution
     */
    public static def make(dist:Dist{rank==1}) : Vector {
        val newVector = new Vector();

        newVector.region       = dist.region;
        newVector.distribution = dist;
        newVector.vec          = Array.make[Double](newVector.distribution);

        return newVector;
    }

    /**
     * Construct a Vector of dimention N, with a custom distribution
     */
    public static def make(newVector:Vector, dist:Dist{rank==1}) : Vector {
        newVector.region       = dist.region;
        newVector.distribution = dist;
        newVector.vec          = Array.make[Double](newVector.distribution);

        return newVector;
    }

 
    /**
     * The size of this matrix
     */
    public global def getSize() : Int = vec.region.max(0)+1;
 
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
        finish foreach((i) in vec.region) vec(i) = 0;
    }

    /**
     * the dot product
     */
    public def dot(b:Vector) : Double {
        var res:Double = 0.0;

        // TODO: this is introduced till a fix for XTENLANG-508 is there
        // finish foreach((i) in vec.region) res += vec(i) * b.vec(i);
        for((i) in vec.region) res += vec(i) * b.vec(i);

        return res;
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val N   = getSize();
        val res = Vector.make(N);

        finish foreach((i) in vec.region) res.vec(i) = vec(i) + b.vec(i);

        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val N   = getSize();
        val res = Vector.make(N);

        finish foreach((i) in vec.region) res.vec(i) = vec(i) - b.vec(i);

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
        
        val n = Vector.make(N);
       
        finish foreach((i) in vec.region) n.vec(i) = vec(i) / mag;
        
        return n;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val N = getSize();
        val n = Vector.make(N);
        
        finish foreach((i) in vec.region) n.vec(i) = -vec(i);
        
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {        
        val N   = getSize();
        val res = Vector.make(N);
        
        finish foreach((i) in vec.region) res.vec(i) = vec(i) * k;
        
        return res;
    }

    public global safe def toString() : String {
         var str : String = "";
         val N = getSize();

         for(var i:Int=0; i<N; i++) {
            str += vec(i) + " ";
         } // end for

         return str;
    }

}

