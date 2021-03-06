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

import x10.io.File;
import x10.io.FileWriter;
import x10.io.FileReader;

/**
 * Evaluate 1E integrals
 *
 * @author: V.Ganesh
 */
public class OneElectronIntegrals(numBasisFunctions:Long) { 
    val basisFunctions:BasisFunctions;
    val hCore:HCore{self.N==numBasisFunctions,self.M==numBasisFunctions};
    val overlap:Overlap{self.N==numBasisFunctions,self.M==numBasisFunctions};
    val roZ:Double;
  
    public def this(bfs:BasisFunctions, mol:Molecule[QMAtom], inpFile:String) {
        property(bfs.getBasisFunctions().size());
        this.basisFunctions = bfs;

        hCore    = new HCore(numBasisFunctions);
        overlap  = new Overlap(numBasisFunctions);

        val jd = JobDefaults.getInstance();
        this.roZ=jd.roZ;
        
        val file = new File(inpFile+".1int");
        var r:FileReader=null;
        try {
            r=new FileReader(file);
            for(i in 0..(hCore.M-1)) {
                for (j in 0..(hCore.N-1)) {
                    overlap(i,j)=r.readDouble();
                    hCore(i,j)=r.readDouble()*roZ;
                }
            }
            Console.OUT.println("1e Ints read from "+inpFile+".1int (Abort the calculation and delete the file if geometry/basis set has changed!)");
        }   catch (ioe:x10.io.IOException) { // include catch (eof:x10.io.EOFException)
            Console.ERR.println(ioe);
            Console.OUT.printf("Calculating 1e Ints.\n");
            compute1E(mol); 
            val w = new x10.io.FileWriter(file);
            for(i in 0..(hCore.M-1)) {
                for (j in 0..(hCore.N-1)) {
                    w.writeDouble(overlap(i,j));
                    w.writeDouble(hCore(i,j)/roZ);
                }
            }
            w.close();
            Console.OUT.println("1e Ints saved to "+inpFile+".1int (Delete the file if geometry/basis set has changed!)");
        } finally {
            if (r!=null) r.close();
        }

    } 

    public def getHCore() = hCore;
    public def getOverlap() = overlap;
    public def getBasisFunctions() = basisFunctions;

    private def compute1E(molecule:Molecule[QMAtom]) : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nat  = molecule.getNumberOfAtoms();
       val atno = new Rail[Double](nat);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0n; i<nat; i++) 
           atno(i) = atms.get(i).species;

       for (i in 0n..(hCore.M-1n)) {
            for (j in 0n..(hCore.N-1n)) {
                val bfi = bfs.get(i);
                val bfj = bfs.get(j);

                val oVal = bfi.overlap(bfj);
                val hVal = bfi.kinetic(bfj)/roZ; /*kinetic scales differently as it is not Coulombic interation*/

                overlap(i,j) = oVal;
                hCore(i,j) = hVal;
                         
                for(var k:Int=0n; k<nat; k++) {
                    val aVal = atno(k) * bfi.nuclear(bfj, atms.get(k).centre);
                    hCore(i,j) += aVal;
                }
            }
       }
    }
}

