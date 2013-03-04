/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 * (C) Copyright Josh Milthorpe 2011-2013.
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

    /**
     * Construct a Vector of dimension N
     */
    public def this(siz:Int) { 
        vec = new Rail[Double](siz);
    }

    /**
     * Construct a Vector from given backing storage
     */
    public def this(vec:Rail[Double]) { 
        this.vec = vec;
    }

    private static def makeUnsafe(size:Int) {
        return new Vector(size);
        // TODO for X10 2.4
        //val store = Unsafe.allocRailUninitialized[Double](size);
        //return new Vector(store);
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

    public @Inline operator this(i0:Int) : Double {
        return vec(i0);
    }

    public @Inline operator this(i0:Int)=(v:Double) : Double {
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
        for(i in 0..(vec.size-1)) {
            res += vec(i) * b.vec(i);
        }
        return res;
    }

    /**
     * add two vectors: this + b 
     */
    public def add(b:Vector) : Vector {
        val res = Vector.makeUnsafe(vec.size);
        for(i in 0..(vec.size-1)) {
            res(i) = vec(i) + b.vec(i);
        }
        return res;
    }

    /**
     * subtract two vectors: this - b
     */
    public def sub(b:Vector) : Vector {
        val res = Vector.makeUnsafe(vec.size);
        for(i in 0..(vec.size-1)) {
            res(i) = vec(i) - b.vec(i);
        }

        return res;
    }

    /**
     * magnitude of this vector
     */
    public def magnitude() : Double {
        var magSquared:Double = 0.0;
        for(i in 0..(vec.size-1)) {
            magSquared += vec(i) * vec(i);
        }
        return Math.sqrt(magSquared);
    }

    /**
     * return a normalized form of this vector
     */
    public def normalize() : Vector {
        val mag = magnitude();
        val res = Vector.makeUnsafe(vec.size);
        for(i in 0..(vec.size-1)) {
            res(i) = vec(i) / mag;
        }
        return res;
    }

    /**
     * return (-1) . V
     */
    public def negate() : Vector {
        val res = Vector.makeUnsafe(vec.size);
        for(i in 0..(vec.size-1)) {
            res(i) = -vec(i);
        }
        return res;
    }

    /**
     * Multiply this vector by a constant k
     */
    public def mul(k:Double) : Vector {
        val res = Vector.makeUnsafe(vec.size);
        for(i in 0..(vec.size-1)) {
            res(i) = vec(i) * k;
        }
        return res;
    }

    /**
     * max norm of this vector
     */
    public def maxNorm() : Double {
        var max:Double = 0.0;
        for(i in 0..(vec.size-1)) {
            max = Math.max(Math.abs(vec(i)), max);
        }
        return max;
    }

    public def toString() : String {
        var str : String = "";

        for(v in vec) {
            str += "" + v + " ";
        }

        return str;
    }
}

