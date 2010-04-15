package x10x.vector;

import x10x.matrix.Matrix;

/**
 * This class represents a N-dimentional Vector
 * (Initial DRAFT)
 * 
 * @author V.Ganesh
 */
public class Vector { 
    global val vec:Array[Double]{rank==1, self.at(this)}; 
    global val region:Region{rank==1, self.at(this)};
    global val distribution:Dist{rank==1, self.at(this)};

    /**
     * Construct a Vector of dimention N, with default block distribution
     */
    public def this(siz:Int) { 
        region       = [0..(siz-1)];
        // distribution = Dist.makeBlock(region, 0) as Dist{rank==1, self.at(this)};
        distribution = Dist.makeConstant(region) as Dist{rank==1, self.at(this)};
        vec          = new Array[Double](region) as Array[Double]{rank==1, self.at(this)};
    }

    /**
     * Construct a Vector from a Matrix
     */
    public def this(mat:Matrix!) {
        this(mat.getRowCount()*mat.getColCount());

        val N = getSize();
         
        var ii:Int = 0;
        val m = mat.getMatrix();

        for((i,j) in m.region)
            vec(ii++) =  m(i, j);
    }

    /**
     * Construct a Vector of dimention N, with a custom distribution
     */
    public def this(dist:Dist{rank==1}) {
        distribution = dist as Dist{rank==1, self.at(this)};
        region       = distribution.region as Region{rank==1, self.at(this)};
        vec          = new Array[Double](region) as Array[Double]{rank==1, self.at(this)};
    }
 
    /**
     * The size of this matrix
     */
    public global def getSize() = vec.region.max(0)+1;
 
    /**
     * Actual data stored
     */
    public def getVector() = vec;

    /**
     * Get associated region
     */ 
    public def region() = region;
 
    /**
     * Get associated distribution
     */
    public def dist() = distribution;

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
        val res = Rail.make[Double](1, (Int)=>0.0);

        // TODO: this is introduced till a fix for XTENLANG-508 is there
        // finish foreach((i) in vec.region) res += vec(i) * b.vec(i);
        finish foreach((i) in vec.region) res(0) += vec(i) * b.vec(i);

        return res(0);
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val N   = getSize();
        val res = new Vector(N);

        finish foreach((i) in vec.region) res.vec(i) = vec(i) + b.vec(i);

        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val N   = getSize();
        val res = new Vector(N);

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
        
        val n = new Vector(N);
       
        finish foreach((i) in vec.region) n.vec(i) = vec(i) / mag;
        
        return n;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val N = getSize();
        val n = new Vector(N);
        
        finish foreach((i) in vec.region) n.vec(i) = -vec(i);
        
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {        
        val N   = getSize();
        val res = new Vector(N);
        
        finish foreach((i) in vec.region) res.vec(i) = vec(i) * k;
        
        return res;
    }

    /**
     * max norm of this vector
     */
    public def maxNorm() : Double {
        val N   = getSize();
        val res = Rail.make[Double](1, (Int)=>0.0);

        finish foreach((i) in vec.region) res(0) = Math.max(Math.abs(vec(i)), res(0));
        
        return res(0);       
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

