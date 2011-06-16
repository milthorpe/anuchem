/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.pme;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;

/**
 * Tests the Distributed Particle Mesh Ewald implementation - basic functions.
 * @author milthorpe
 */
public class TestPMEFunctions extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public def this(numAtoms : Int) {
        super(numAtoms);
    }

    public static def main(args : Array[String](1)) {
        new TestPMEFunctions(2).test();
    }

    public def test() {
        val atoms = generateAtoms();
        val edges = [Vector3d(SIZE, 0.0, 0.0), Vector3d(0.0, SIZE, 0.0), Vector3d(0.0, 0.0, SIZE)];
        val gridSize = new Array[Int](3, 4);
        val pme = new PME(edges, gridSize, atoms, 4, 0.35, 9.0);
        val v = Vector3d(0.6, -0.4, 0.33);
        val r = pme.getScaledFractionalCoordinates(v);
        Console.OUT.println(r);

        val m404 = pme.bSpline(4, 0);
        Console.OUT.println("M_4(3) = " + m404);

/*
        val Q = pme.getGriddedCharges();
        var nonZero : Int = 0;
        for (p in Q) {
            val q = at(Q.dist(p)) {Q(p)};
            if (q != Complex.ZERO) {
                nonZero++;
                //Console.OUT.println(p + " = " + q);
            }
        }
        Console.OUT.println("nonZero = " + nonZero);
*/

        val B = pme.getBArray();
        var sum : Double = 0.0;
        for (p in B) {
            val b = at(B.dist(p)) {B(p)};
            sum += B(p);
            if (b != 0.0) {
                Console.OUT.println("B" + p + " = " + b);
            }
        }
        Console.OUT.println("sum B = " + sum);

        val C = pme.getCArray();
        for (p in C) {
            val c = at(C.dist(p)) {C(p)};
            if (c != 0.0) {
                Console.OUT.println("C" + p + " = " + c);
            }
        }

        val BdotC = DistArray.make[Double](B.dist);
	    finish ateach(p in B.dist) {
		    BdotC(p) = B(p) * C(p);
	    }
        for (p in BdotC) {
            at (BdotC.dist(p)) {
                if (BdotC(p) != 0.0) {
                    Console.OUT.println("BdotC" + p + " = " + BdotC(p));
                }
            }
        }
    }
}

