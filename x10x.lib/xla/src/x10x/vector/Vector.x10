/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 * (C) Copyright Josh Milthorpe 2011.
 */
package x10x.vector;

import x10.compiler.Inline;
import x10x.matrix.Matrix;

/**
 * This class represents a N-dimensional vector
 * TODO distributed Vector
 * 
 * @author V.Ganesh, milthorpe
 */
public class Vector { 
    val vec:Rail[Double];
    property region() = vec.region;

    /**
     * Construct a Vector of dimention N
     */
    public def this(siz:Int) { 
        vec = new Array[Double](siz);
    }

    /**
     * Construct a Vector of dimention N
     */
    public def this(region:Region(1){rect,zeroBased,rail}) { 
        vec = new Array[Double](region);
    }

    /**
     * Construct a Vector row-packed from a Matrix
     */
    public def this(mat:Matrix) {
        this(mat.getRowCount()*mat.getColCount());
         
        var ii:Int = 0;
        val m = mat.getMatrix();

        for([i,j] in m.region)
            vec(ii++) =  m(i, j);
    }
 
    /**
     * The size of this vector
     */
    public def getSize() = vec.size;
 
    /**
     * Actual data stored
     */
    public def getVector() = vec;

    public @Inline operator this(i0:int) : Double {
        return vec(i0);
    }

    public @Inline operator this(i0:int)=(v:Double) : Double {
        vec(i0) = v;
        return v;
    }

    /**
     * make this Vector a null vector
     */
    public def makeZero() : void {
        vec.clear();
    }

    /**
     * the dot product
     */
    public def dot(b:Vector) : Double {
        var res : Double = 0.0;
        for([i] in vec) {
            res += vec(i) * b.vec(i);
        }
        return res;
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val res = new Vector(this.region);

        val plus = ((a : Double, b : Double) => a + b);
        vec.map[Double,Double](res.vec, b.vec, plus);

        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val res = new Vector(this.region);
        val minus = ((a : Double, b : Double) => a - b);
        vec.map[Double,Double](res.vec, b.vec, minus);

        return res;
    }

    /**
     * magnitude of this vector
     */
    public def magnitude() : Double {
        val magSquared = ((res : Double, a : Double) => res + a*a);
        return Math.sqrt(vec.reduce(magSquared, 0.0));
    }

    /**
     * return a normalized form of this vector
     */
    public def normalize() : Vector {
        val mag = magnitude();
        val res = new Vector(this.region);
        val norm = ((a : Double) => a / mag);
        vec.map[Double](res.vec, norm);
        return res;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val n = new Vector(this.region);
        val neg = ((a : Double) => -a);
        vec.map[Double](n.vec, neg);
        return n;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {
        val res = new Vector(this.region);
        val mulK = ((a : Double) => a * k);
        vec.map[Double](res.vec, mulK);
        return res;
    }

    /**
     * max norm of this vector
     */
    public def maxNorm() : Double {
        val max = ((res : Double, a : Double) => Math.max(Math.abs(a), res));
        return vec.reduce(max, -1.0); 
    }

    public def toString() : String {
        var str : String = "";

        for([i] in vec.region) {
            str += "" + vec(i) + " ";
        }

        return str;
    }
}

