/**
 * OneElectronIntegrals.x10
 *
 * Evaluate 1E integrals
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

public class OneElectronIntegrals { 
    var basisFunctions:BasisFunctions;
    var molecule:Molecule;
    var hCore:HCore;
    var overlap:Overlap;

    public def this(bfs:BasisFunctions, mol:Molecule) { 
       this.basisFunctions = bfs;
       this.molecule = mol;
 
       compute1E();
    } 

    public def getHCore() : HCore = hCore;
    public def getOverlap() : Overlap = overlap;
    public def getBasisFunctions() : BasisFunctions = basisFunctions;
    public def getMolecule() : Molecule = molecule;

    protected def compute1E() : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nbf  = bfs.size();
       hCore    = new HCore(nbf);
       overlap  = new Overlap(nbf);
       val nat  = molecule.getNumberOfAtoms();
       val atno = Array.make[Double]([0..nat]);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0; i<nat; i++) atno(i) = ai.getAtomicNumber(atms.get(i));

       val ovr = overlap.getMatrix();
       val h   = hCore.getMatrix();

       // TODO : x10 - parallel
       for(var i:Int=0; i<nbf; i++) {
          var bfi:ContractedGaussian = bfs.get(i);

          for(var j:Int=0; j<nbf; j++) {
             var bfj:ContractedGaussian = bfs.get(j);
             
             ovr(i,j) = bfi.overlap(bfj); 
             h(i,j)   = bfi.kinetic(bfj);
 
             for(var k:Int=0; k<nat; k++)
                h(i,j) += atno(k) * bfi.nuclear(bfj, atms.get(k));
          } // end for
       } // end for
    }
}

