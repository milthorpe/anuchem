/**
 * The Gaussian elemination solver for A x = B
 *
 * Node: Mostly lifted from MeTA Studio code, except no Matrix singularity 
 * and upper - triangular checks made!
 * 
 * @author  V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */

package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

public class GaussianElimination extends LinearEquationSolver {

    private var row:Array[Int]{rank==1, self.at(this)};
    private var a:Array[Double]{rank==2, self.at(this)};
    private var x:Array[Double]{rank==1, self.at(this)};
    
    public def this() {
    }

    private var n:Int, n1:Int;

    public def findSolution(matrixA:Matrix!, vectorB:Vector!) : Vector! {
        val N = matrixA.getRowCount();
        val m = new Matrix(Dist.make([0..N-1, 0..N])) as Matrix!;
        a = m.getMatrix() as Array[Double]{rank==2, self.at(this)};
        x = vectorB.getVector() as Array[Double]{rank==1, self.at(this)};

        var i:Int, j:Int;

        val ta = matrixA.getMatrix();
        for(i=0; i<N; i++)
            for(j=0; j<N; j++)
                a(i,j) = ta(i,j);

        for(i=0; i<N; i++)
            a(i,N) = x(i);        
        
        n  = N-2;
        n1 = N-1;
        
        row = Array.make[Int]([0..n]) as Array[Int]{rank==1, self.at(this)};
        
        // initilize row vector
        for(j=1; j<=n; j++)
            row(j) = j;

        // first make the the matrix upper triangular
        upperTriangularize(true);

        // now use back substitution to get the solution
        // first solution from back :)
        x(n) = a(row(n),n1) / a(row(n),n);

        var sum:Double = 0.0;
        var k:Int, c:Int;
        for(k=(n-1); k>=0; k--) {
            sum = 0.0;
            for(c=(k+1); c<=n; c++) {
                sum += (a(row(k),c) * x(c)); // compute next solutions
            } // end of inner loop

            // the actual next solution
            x(k) = (a(row(k),n1) - sum) / a(row(k),k);
        } // end outer for

        return vectorB;
    }

    /** Upper triangularize the matrix */
    private def upperTriangularize(doOneScale:Boolean) : Boolean {
        // check if matrix is already upper triangular
        // this check may be removed in future
        /** 
         // TODO: skipping this!
        if (matrixA.isUpperTriangular()) {
            return;
        } // end if
        **/

        for (var p:Int=0; p<=n1; p++) {
            simplePivot(p); // apply simple pivoting
            oneScale(p, doOneScale);    // apply simple 1-scaling
        } // end for

        /** 
         // TODO: Should ideally be there ...
        if (matrixA.isSingular(n, row)) {
            return false;
        } // end if
        **/

        return true;
    }

    /**
     * simplePivot() - This method does simple pivoting if the a[j][j]th element
     *                 is zero. It first finds the a[i][j]th element which is
     *                 numerically greater than a[j][j] and i>j and then
     *                 interchanges the ith and jth rows.
     * @param p - The pth iteration in Gaussian elemination.
     */
    public def simplePivot(p:Int) : Boolean {
        var temp:Int = 0;

        if (p >= n) {
            return true;
        } // end if

        // if concerned row's first element is unity
        // or more do not pivot
        if (Math.abs(a(row(p),p)) >= 1) return true;

        var k:Int;
        for(k=(p+1); k<=n; k++) {
            // check for keeping things near unity
            if ((Math.abs(a(row(k),p))-1) < (Math.abs(a(row(p),p))-1)) {
                // switch the indices to represent a
                // row interchange so as to avoide the overhead
                // of actually moving the row elements  :)
                temp   = row(p);
                row(p) = row(k);
                row(k) = temp;
            } // end if
        } // end for

        /** 
         // TODO: Should ideally be there ...
         // check if singular
        if (matrixA.isSingular(p, row)) {
           return false;
        } // end if
        */

        return true;
    } // end of method simplePivot()

    /**
     * oneScale() - This method makes the a[j][j]th entry at the jth iteration
     *              close to unity by deviding each element of the jth row
     *              by a[j][j].
     * @param p - The jth iteration in Gaussian elemination.
     *        boolean scale - Scale to one or not
     */
    private def oneScale(p:Int, scale:Boolean) {
        var m:Double = 0.0;

        if (p >= n) {
            return;
        } // end if
 
        var k:Int, c:Int;

        for(k=(p+1); k<=n; k++) {
            // compute the multiplicity factor
            m = a(row(k),p) / a(row(p),p);

            for(c=(p+1); c<=n1; c++) {
                a(row(k),c) -= (m * a(row(p),c));
            } // end inner for

            a(row(k),p) = m;
        } // end outer for

        if (scale) {
            a(row(p),p) = 1.0; // diagonal element is now 1
        } // end if
        
        for(c=(p+1); c<=n; c++) {
            a(row(c),p) = 0.0; // all entries below diagonal are 0
        } // end for
    } // end of method oneScale()
}


