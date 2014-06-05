/*  This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2014.
 */
package au.edu.anu.util;

import x10.compiler.Ifdef;
import x10.compiler.Native;
import x10.compiler.NativeCPPInclude;
import x10.io.IOException;
import x10.regionarray.Dist;

/** 
 * Provides methods to query and modify the execution environment ANUChem.
 */
@NativeCPPInclude("mkl_math.h")
@NativeCPPInclude("omp.h")
public class ExecutionEnvironment {
    @Native("c++", "omp_get_num_threads()") private native static def ompGetNumThreads():Int;
    @Native("c++", "omp_set_num_threads(#a)") private native static def ompSetNumThreads(a:Int):void;
    @Native("c++", "mkl_get_max_threads()") private native static def mklGetMaxThreads():Int;
    @Native("c++", "MKL_Set_Num_Threads(#a)") private native static def mklSetNumThreads(a:Int):void;

    /** Change the number of threads used for BLAS and LAPACK calls at this place. */
    public static def setBlasThreads(numThreads:Int) {
        @Ifdef("__MKL__") { // not working?  better use -genv OMP_NUM_THREADS 4
            val t1=mklGetMaxThreads();
            val o1=ompGetNumThreads();
            ompSetNumThreads(numThreads);
            mklSetNumThreads(numThreads);
            val t2=mklGetMaxThreads();
            val o2=ompGetNumThreads();
            @Ifdef("__DEBUG__") { Console.OUT.println(here + ", mklGetMaxThreads() was " + t1 + " and is now set to " + t2 + " thread(s)."
                            + " ompGetNumThreads() was " + o1 + " and is now set to " + o2 + " thread(s)."); }   
        }
    }

    /** Print threading environment variables at each place. */
    public static def printThreadingVariables() {
        finish ateach(place in Dist.makeUnique()) {
            val hostname=Runtime.execForRead("uname -n").readLine();
            // should be System.getenv("X10_NPLACES") - XTENLANG-256
            val np = System.getenv().getOrElse("X10_NPLACES", "");
            val nt = System.getenv().getOrElse("X10_NTHREADS", "");
            val gt = System.getenv().getOrElse("GOTO_NUM_THREADS", "");
            val omp = System.getenv().getOrElse("OMP_NUM_THREADS", "");
            Console.OUT.println(here + ", Runtime.NTHREADS=" + Runtime.NTHREADS + ", uname -n=" + hostname + ", X10_NPLACES="+np+", X10_NTHREADS="+nt+", GOTO_NUM_THREADS="+gt+ ", OMP_NUM_THREADS="+omp );
            Console.OUT.flush();
        }
   }
}
