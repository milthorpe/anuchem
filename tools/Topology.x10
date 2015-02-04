/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2013.
 */

import x10.xrx.Runtime;

/**
 * This class prints the hostname and number of threads for each place in
 * the X10 system.  It is useful for understanding the topology of the system.
 */
public class Topology {
    public static def main(args:Rail[String]){
        for (place in Place.places()) at(place) {
            val hostnameReader = Runtime.execForRead("uname -n"); 
            val hostname = hostnameReader.readLine();
            Console.OUT.println(here + " executing on " + Runtime.getName() + " with " + Runtime.NTHREADS + " threads");

/*
            // useful for debugging X10 environment

            val npReader  = Runtime.execForRead("echo $X10_NPLACES");
            val np = npReader.readLine();
            Console.OUT.println("X10_NPLACES="+np);

            val ntReader  = Runtime.execForRead("echo $X10_NTHREADS");
            val nt = ntReader.readLine();
            Console.OUT.println("X10_NTHREADS="+nt);

            // useful for BLAS libraries (Intel-MKL, GotoBLAS2)

            val gtReader  = Runtime.execForRead("echo $GOTO_NUM_THREADS");
            val gt = gtReader.readLine();
            Console.OUT.println("GOTO_NUM_THREADS="+gt);

            val ompReader  = Runtime.execForRead("echo $OMP_NUM_THREADS");
            val omp = ompReader.readLine();
            Console.OUT.println("OMP_NUM_THREADS="+omp);
*/
            
        }
    }
}
