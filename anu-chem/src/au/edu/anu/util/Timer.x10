/** lifted from HPC benchmarks example */

package au.edu.anu.util;

final public class Timer {
    public val total:Array[Long](1);
    public val count:Array[Long](1);

    public def this(n:Int) {
        total = new Array[Long](n);
        count = new Array[Long](n);
    }

    public def start(id:Int) { total(id) -= System.nanoTime(); }
    public def clear(id:Int) { total(id) = 0; }
    public def stop(id:Int) { total(id) += System.nanoTime(); count(id)++; }
    public def mean(id:Int) = total(id) / count(id);
}

