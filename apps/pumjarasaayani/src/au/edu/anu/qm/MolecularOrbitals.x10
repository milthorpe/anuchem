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

import x10.matrix.DenseMatrix;

/**
 * Represents MOs in a HF-SCF
 *
 * @author: V.Ganesh
 */
public class MolecularOrbitals extends DenseMatrix {
    private var orbitalEnergies:Rail[Double];

    public def this(n:Int) {
        super(n, n);
    }

    public def getOrbitalEnergies() : Rail[Double] = orbitalEnergies;
 
    public def compute(theMat:DenseMatrix(this.N,this.N), overlap:Overlap{self.M==this.M,self.N==this.N}) : void {
        val x = overlap.getSHalf();
        val a = new DenseMatrix(x.M, theMat.N);
        a.multTrans((x % theMat), x); // a = x.theMat.x^T

        val diag = new GMLDiagonalizer();

        diag.diagonalize(a);
        orbitalEnergies = diag.getEigenValues().d;

        super.mult(diag.getEigenVectors(), x);
    }
}

