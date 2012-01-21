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
package au.edu.anu.qm;
import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

/**
 * Evaluate 1E integrals
 *
 * @author: V.Ganesh
 */
public class OneElectronIntegrals(numBasisFunctions:Int) { 
    val basisFunctions:BasisFunctions;
    val hCore:HCore{self.N==numBasisFunctions,self.M==numBasisFunctions};
    val overlap:Overlap{self.N==numBasisFunctions,self.M==numBasisFunctions};
    val roZ:Double;
  
    public def this(bfs:BasisFunctions, mol:Molecule[QMAtom]) {
        property(bfs.getBasisFunctions().size());
        this.basisFunctions = bfs;

        hCore    = new HCore(numBasisFunctions);
        overlap  = new Overlap(numBasisFunctions);

        val jd = JobDefaults.getInstance();
        this.roZ=jd.roZ;

        compute1E(mol);
    } 

    public def getHCore() = hCore;
    public def getOverlap() = overlap;
    public def getBasisFunctions() = basisFunctions;

    private def compute1E(molecule:Molecule[QMAtom]) : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nat  = molecule.getNumberOfAtoms();
       val atno = new Array[Double](nat);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0; i<nat; i++) 
           atno(i) = ai.getAtomicNumber(atms.get(i));

       for([i, j] in 0..(hCore.M-1)*0..(hCore.N-1)) {
              val bfi = bfs.get(i);
              val bfj = bfs.get(j);

              val oVal = bfi.overlap(bfj);
              val hVal = bfi.kinetic(bfj)/roZ; /*kinetic scales differently as it is not Coulombic interation*/

              overlap(i,j) = oVal;
              hCore(i,j) = hVal;
                         
              for(var k:Int=0; k<nat; k++) {
                  val aVal = atno(k) * bfi.nuclear(bfj, atms.get(k).centre);
                  hCore(i,j) += aVal;
              }
           
       }
    }
}

