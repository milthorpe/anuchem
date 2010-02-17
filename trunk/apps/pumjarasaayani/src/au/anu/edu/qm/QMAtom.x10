/**
 * QMAtom.x10
 *
 * This class represents an atom for the purposes of 
 * Quantum Chemistry simulations
 *
 * @author: V.Ganesh
 */
package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.vector.Point3d;
import au.edu.anu.chem.Atom;

public class QMAtom extends Atom {
    public def this(symbol : String, centre : Point3d!) {
        super(symbol, centre);
    }

    public def this(centre : Point3d!) {
        super(centre);
    }

    global var basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;
}

