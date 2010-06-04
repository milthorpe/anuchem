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

    val dummy:Boolean;

    public def this(symbol : String, centre : Point3d) {
        super(symbol, centre);
        index = 0;
        dummy = false;        
    }

    public def this(centre : Point3d) {
        super(centre);
        index = 0;
        dummy = false;
    }

    public def this(symbol : String, centre : Point3d, dummy : Boolean) {
        super(symbol, centre);
        index = 0;
        this.dummy = dummy;
    }

    public def isDummy() = dummy;
    
    var basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;
}

