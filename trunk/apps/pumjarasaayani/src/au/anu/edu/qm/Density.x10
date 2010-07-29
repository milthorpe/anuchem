/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.matrix.Matrix;

/**
 * Density.x10
 *
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Density extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    public def compute(noOfOccupancies:Int, mos:MolecularOrbitals!) : void {
        // unmute();

        // construct it from the MOs .. C*C'
        val N = mos.getRowCount();
        val dVector:Matrix! = new Matrix(noOfOccupancies, N) as Matrix!;

        val dMat = dVector.getMatrix();
        val mosMat = mos.getMatrix();
        for(var i:Int=0; i<noOfOccupancies; i++)
            for(var j:Int=0; j<N; j++) 
              dMat(i, j) = mosMat(i, j);

        val res = dVector.transpose().mul(dVector).getMatrix();

        val thisMat = getMatrix();
        for(val(i, j) in thisMat.region)
           thisMat(i, j) = res(i, j);

        // mute();
    }
}

