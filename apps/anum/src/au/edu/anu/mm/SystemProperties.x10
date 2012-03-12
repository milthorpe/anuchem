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
    public val numAtoms:Int;
    public val oneParticleFunctions:Rail[OneParticleFunction];
    
    public def this(numAtoms:Int, oneParticleFunctions:Rail[OneParticleFunction]) {
        this.numAtoms = numAtoms;
        this.oneParticleFunctions = oneParticleFunctions;
    }

    public def calculatePropertySums(atoms: DistArray[Rail[MMAtom]](1)):Rail[Pair[String,Double]] {
        val totals = new Accumulator[Rail[Double]](RailSumReducer(oneParticleFunctions.size));
        finish ateach(place in atoms) {
            val atomsHere = atoms(place);
            val totalHere = calculatePropertySums(atomsHere);
            totals <- totalHere;
        }
        
        val raw = totals();
        val results = new Rail[Pair[String,Double]](raw.size);
        for (i in 0..(results.size-1)) {
            results(i) = Pair[String,Double](oneParticleFunctions(i).first, raw(i) / numAtoms);
        }
        return results;
    }

    public def calculatePropertySums(atoms:Rail[MMAtom]):Rail[Double] {
        val totals = new Rail[Double](oneParticleFunctions.size);
        for ([p] in atoms) {
            for ([i] in oneParticleFunctions) {
                totals(i) += oneParticleFunctions(i).second(atoms(p));
            }
        }
        return totals;
    }

    static struct RailSumReducer(size:Int) implements Reducible[Rail[Double]] {
        public def this(size:Int) { property(size); }
        public def zero() = new Rail[Double](size);
        public operator this(a:Rail[Double], b:Rail[Double]) {
            val result = new Rail[Double](b.size);
            for (i in 0..(b.size-1)) {
                result(i) = a(i) + b(i);
            }
            return result;
        }
    }
}

public static type OneParticleFunction = Pair[String,(a:MMAtom) => Double];
public static type TwoParticleFunction = Pair[String,(a:MMAtom,b:MMAtom) => Double];

