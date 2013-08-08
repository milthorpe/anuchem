/** lifted from HPC benchmarks example */

package au.edu.anu.util;

final public class Timer {
    public val total:Rail[Long];
    public val count:Rail[Long];

    public def this(n:Long) {
        total = new Rail[Long](n);
        count = new Rail[Long](n);
    }

    public def clear() { total.clear(); count.clear(); }

    public def start(id:Long) { total(id) -= System.nanoTime(); }
    public def clear(id:Long) { total(id) = 0L; count(id) = 0L; }
    public def stop(id:Long) { total(id) += System.nanoTime(); count(id)++; }
    public def mean(id:Long) {
        if (count(id) == 0L) return 0L;
        else return total(id) / count(id);
    }
}

