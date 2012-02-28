/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.util.Pair;
import au.edu.anu.chem.mm.MMAtom;

/**
 * This class manages calculation of dynamical and ensemble system properties.
 */
public class SystemProperties {
    public val oneParticleFunctions:Rail[OneParticleFunction];
    
    public def this(oneParticleFunctions:Rail[OneParticleFunction]) {
        this.oneParticleFunctions = oneParticleFunctions;
    }

    public def calculateExpectationValues(atoms: DistArray[Rail[MMAtom]](1)):Rail[Pair[String,Double]] {
        val totals = new Accumulator[Rail[Pair[Int,Double]]](RailSumReducer(oneParticleFunctions.size));
        finish ateach(place in atoms) {
            val totalHere = new Rail[Pair[Int,Double]](oneParticleFunctions.size);
            val atomsHere = atoms(place);
            for ([p] in atomsHere) {
                for ([i] in oneParticleFunctions) {
                    val atomVal = oneParticleFunctions(i).second(atomsHere(p));
                    val currentVal = totalHere(i);
                    totalHere(i) = Pair[Int,Double](currentVal.first+1, currentVal.first+atomVal);
                }
            }
            totals <- totalHere;
        }
        
        val raw = totals();
        val results = new Rail[Pair[String,Double]](raw.size);
        for (i in 0..(results.size-1)) {
            results(i) = Pair[String,Double](oneParticleFunctions(i).first, raw(i).second / raw(i).first);
        }
        return results;
    }

    static struct RailSumReducer(size:Int) implements Reducible[Rail[Pair[Int,Double]]] {
        public def this(size:Int) { property(size); }
        public def zero() = new Rail[Pair[Int, Double]](size);
        public operator this(a:Rail[Pair[Int, Double]], b:Rail[Pair[Int, Double]]) {
            val result = new Rail[Pair[Int, Double]](b.size);
            for (i in 0..(b.size-1)) {
                result(i) = Pair[Int,Double](a(i).first + b(i).first, a(i).second + b(i).second);
            }
            return result;
        }
    }
}

public static type OneParticleFunction = Pair[String,(a:MMAtom) => Double];
public static type TwoParticleFunction = Pair[String,(a:MMAtom,b:MMAtom) => Double];

