/**
 * OneElectronIntegrals.x10
 *
 * Evaluate 1E integrals
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;
import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

public class OneElectronIntegrals { 
    global val basisFunctions:BasisFunctions{self.at(this)};
    global val hCore:HCore{self.at(this)};
    global val overlap:Overlap{self.at(this)};

    public def this(bfs:BasisFunctions!, mol:Molecule[QMAtom]!) { 
       this.basisFunctions = bfs;

       val nbf  = basisFunctions.getBasisFunctions().size();
       hCore    = HCore.make(new HCore(), nbf) as HCore{self.at(this)};
       overlap  = Overlap.make(new Overlap(), nbf) as Overlap{self.at(this)};

       compute1E(mol);
    } 

    public def getHCore() = hCore;
    public def getOverlap() = overlap;
    public def getBasisFunctions() = basisFunctions;

    protected def compute1E(molecule:Molecule[QMAtom]!) : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nbf  = bfs.size();
       val nat  = molecule.getNumberOfAtoms();
       val atno = Array.make[Double]([0..nat]);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0; i<nat; i++) atno(i) = ai.getAtomicNumber(atms.get(i));

       val ovr = overlap.getMatrix();
       val h   = hCore.getMatrix();

       // TODO : x10 - parallel
       for(var i:Int=0; i<nbf; i++) {
          var bfi:ContractedGaussian{self.at(this)} = bfs.get(i);

          for(var j:Int=0; j<nbf; j++) {
             var bfj:ContractedGaussian{self.at(this)} = bfs.get(j);
             
             ovr(i,j) = bfi.overlap(bfj); 
             h(i,j)   = bfi.kinetic(bfj);

             // x10.io.Console.OUT.println("K = " + i + ", " + j + " = " + h(i, j));

             for(var k:Int=0; k<nat; k++)
                h(i,j) += atno(k) * bfi.nuclear(bfj, atms.get(k).centre);
          } // end for
       } // end for
    }
}

