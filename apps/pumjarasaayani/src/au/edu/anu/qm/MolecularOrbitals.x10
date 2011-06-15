/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10x.matrix.Matrix;
import x10x.xla.JacobiDiagonalizer;

/**
 * Represents MOs in a HF-SCF
 *
 * @author: V.Ganesh
 */
public class MolecularOrbitals extends Matrix {
    private var orbitalEnergies:Rail[Double];

    public def this(n:Int) {
        super(n);
    }

    public def getOrbitalEnergies() : Rail[Double] = orbitalEnergies;
 
    public def compute(theMat:Matrix, overlap:Overlap) : void {
        val x = overlap.getSHalf();
        val a = theMat.similarityTransform(x);
        val diag = new NativeDiagonalizer();

        diag.diagonalize(a);
        orbitalEnergies = diag.getEigenValues().getVector();

        this.mulInPlace(diag.getEigenVectors(), x);
    }
}

