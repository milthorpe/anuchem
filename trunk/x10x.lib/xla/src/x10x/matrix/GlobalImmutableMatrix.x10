package x10x.matrix;

import x10x.xla.*;
import x10x.vector.Vector;

import x10.array.Array;

/**
 * prototype, programmer controlled global "immutable" array for HF
 *
 * usage:
 *
 *      val ga = new GlobalImmutableMatrix(2);
 *      ga.make();
 *
 *      ga.set(0,0,1);
 *      ga.set(1,0,2);
 *      ga.set(0,1,3);
 *      ga.set(1,1,4);
 *
 *      ga.mute();
 *
 *      for(place in Place.places) {
 *         at(place) {
 *              Console.OUT.println(place + ": " +  ga.get(0,0));
 *              Console.OUT.println(place + ": " +  ga.get(1,0));
 *              Console.OUT.println(place + ": " +  ga.get(0,1));
 *              Console.OUT.println(place + ": " +  ga.get(1,1));
 *         }
 *      }
 *
 *      ga.unmute();
 *
 *      ga.set(0,0,4);
 *      ga.set(1,0,3);
 *      ga.set(0,1,2);
 *      ga.set(1,1,1);
 *
 *      ga.mute();
 *
 *      for(place in Place.places) {
 *         at(place) {
 *              val ar = ga.get();
 *              Console.OUT.println(place + ": " +  ar(0,0));
 *              Console.OUT.println(place + ": " +  ar(1,0));
 *              Console.OUT.println(place + ": " +  ar(0,1));
 *              Console.OUT.println(place + ": " +  ar(1,1));
 *         }
 *      }
 *
 *
 *
 * @author V. Ganesh
 */

public class GlobalImmutableMatrix extends Matrix {

    global val replicatedArray:DistArray[Array[Double]{rect,rank==2}](1);

    global val size:Int;

    global val muted:Rail[Boolean];

    /** Create a new instance of matrix which is square */
    public def this(size:Int) {
        super(size);
        this.size = size;
        	
        replicatedArray = DistArray.make[Array[Double]{rect,rank==2}](Dist.makeUnique());        

        muted = Rail.make[Boolean](1);
    }
 
    public def this(row:Int, col:Int) {
        // TODO: this is unsupported
        super(row, col);

        this.size = 0;
        replicatedArray = DistArray.make[Array[Double]{rect,rank==2}](Dist.makeUnique());

        muted = Rail.make[Boolean](1);
    }

    public def this(dist:Dist{rank==2}) {
        // TODO: this is unsupported
        super(dist);

        this.size = 0;
        replicatedArray = DistArray.make[Array[Double]{rect,rank==2}](Dist.makeUnique());

        muted = Rail.make[Boolean](1);
    }

    /** initialise memory */
    public def make() {                
        finish ateach(idx in replicatedArray.dist)
           replicatedArray(idx) = new Array[Double]([0..size-1,0..size-1]);           
    }

    /**
     * set a value in the matrix, sets are always local. so if muted is 
     * true, then this method simply returns false
     */
    public global def set(i:Int, j:Int, val:Double) : Boolean {        
        if (muted(0)) return false;

        mat(i,j) = val;
        return true;
    }

    /**
     * set all the values in the matrix, sets are always local. so if muted is 
     * true, then this method simply returns false
     */
    public global def set(ar:Array[Double]{rank==2}) : Boolean {
        if (muted(0)) return false;

        ar.copyFrom(mat);
        
        return true;
    }

    /**
     * get a value from the matrix.
     */ 
    public global def get(i:Int, j:Int) : Double {        
        val isMuted = at(muted) { return muted(0); };        
        if (!isMuted) { return mat(i,j); }
        else          { return replicatedArray(here.id)(i,j); }
    }

    /**
     * get all the values of the matrix
     */
    public global def get() : Array[Double]{rank==2} {
        val isMuted = at(muted) { return muted(0); };        
        if (!isMuted) { return mat; }
        else          { return replicatedArray(here.id); }
    }

    /** 
     * mute this matrix.
     * two things happen: first, the values of current array are replicated
     * for all places. second, any subsequent calls to set() are ignored
     * untill unmute() is called.
     */
    public def mute() {       
       muted(0) = true;

       // TODO :
       // 1. possibility of using MPI collectives for broadcast (under the hood)
       // 2. This is a bit weird implementation. This is actually not a
       //    broadcast, but each place "reading" from Place 0
        val myMat = mat;
        finish ateach(idx in replicatedArray.dist) {
            val repArr = replicatedArray(idx);
            repArr.copyFrom(myMat);
        }       
    }

    /** 
     * unmute the matrix, allowing for set() to work.
     */
    public def unmute() {
       muted(0) = false;
    }
}

