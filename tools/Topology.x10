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

/**
 * This class prints the hostname and number of threads for each place in
 * the X10 system.  It is useful for understanding the topology of the system.
 */
public class Topology {
    public static def main(args:Rail[String]){
        for(place in PlaceGroup.WORLD) at(place) {
            val hostnameReader = Runtime.execForRead("uname -n"); 
            val hostname = hostnameReader.readLine();
            Console.OUT.println(here + " executing on " + hostname + " with " + Runtime.NTHREADS + " threads");
        }
    }
}
