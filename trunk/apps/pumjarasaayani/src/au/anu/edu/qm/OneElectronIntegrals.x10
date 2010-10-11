/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;
import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

/**
 * OneElectronIntegrals.x10
 *
 * Evaluate 1E integrals
 *
 * @author: V.Ganesh
 */
public class OneElectronIntegrals { 
    val basisFunctions:BasisFunctions;
    val hCore:HCore;
    val overlap:Overlap;

    public def this(bfs:BasisFunctions, mol:Molecule[QMAtom]) { 
       this.basisFunctions = bfs;

       val nbf  = basisFunctions.getBasisFunctions().size();
       hCore    = new HCore(nbf);
       overlap  = new Overlap(nbf);

       compute1E(mol);
    } 

    public def getHCore() = hCore;
    public def getOverlap() = overlap;
    public def getBasisFunctions() = basisFunctions;

    private def compute1E(molecule:Molecule[QMAtom]) : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nbf  = bfs.size();
       val nat  = molecule.getNumberOfAtoms();
       val atno = new Array[Double](nat);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0; i<nat; i++) 
           atno(i) = ai.getAtomicNumber(atms.get(i));

       val ovr = overlap.getMatrix();
       val h   = hCore.getMatrix();

       finish for([i, j] in h) {
              val bfi = bfs.get(i);
              val bfj = bfs.get(j);

              val oVal = bfi.overlap(bfj);
              val hVal = bfi.kinetic(bfj);

              ovr(i,j) = oVal; 
                h(i,j) = hVal; 

              for(var k:Int=0; k<nat; k++) {
                  val aVal = atno(k) * bfi.nuclear(bfj, atms.get(k).centre);

                  h(i,j) += aVal; 
              }
       }
    }
}

