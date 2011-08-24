/** lifted from HPC benchmarks example */

package au.edu.anu.util;

final public class Timer {
    public val total:Rail[Long];
    public val count:Rail[Long];

    public def this(n:Int) {
        total = new Rail[Long](n);
        count = new Rail[Long](n);
    }

    public def start(id:Int) { total(id) -= System.nanoTime(); }
    public def clear(id:Int) { total(id) = 0; count(id) = 0; }
    public def stop(id:Int) { total(id) += System.nanoTime(); count(id)++; }
    public def mean(id:Int) = total(id) / count(id);
}

