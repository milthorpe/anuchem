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
    global val basisFunctions:BasisFunctions!;
    global val hCore:HCore!;
    global val overlap:Overlap!;

    public def this(bfs:BasisFunctions!, mol:Molecule[QMAtom]!) { 
       this.basisFunctions = bfs;

       val nbf  = basisFunctions.getBasisFunctions().size();
       hCore    = new HCore(nbf) as HCore!;
       overlap  = new Overlap(nbf) as Overlap!;

       compute1E(mol);
    } 

    public def getHCore() = hCore;
    public def getOverlap() = overlap;
    public def getBasisFunctions() = basisFunctions;

    protected def compute1E(molecule:Molecule[QMAtom]!) : void {
       val bfs  = basisFunctions.getBasisFunctions();
       val nbf  = bfs.size();
       val nat  = molecule.getNumberOfAtoms();
       val atno = Rail.make[Double](nat);
       val ai   = AtomInfo.getInstance();
       val atms = molecule.getAtoms();

       for(var i:Int=0; i<nat; i++) 
           atno(i) = ai.getAtomicNumber(atms.get(i));

       val ovr = overlap.getMatrix();
       val h   = hCore.getMatrix();

       finish foreach(plc in h.dist.places()) {
             for(val(i, j) in h.dist.get(plc)) {
                 val bfi = bfs.get(i);
                 val bfj = bfs.get(j);

                 val oVal = bfi.overlap(bfj);
                 val hVal = bfi.kinetic(bfj);

                 at(plc) { ovr(i,j) = oVal; h(i, j) = hVal; }

                 for(var k:Int=0; k<nat; k++) {
                    val aVal = atno(k) * bfi.nuclear(bfj, atms.get(k).centre);

                    at(plc) { h(i,j) += aVal; }
                 }
             }
       }
    }
}

