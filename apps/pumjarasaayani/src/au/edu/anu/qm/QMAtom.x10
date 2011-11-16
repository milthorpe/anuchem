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

import x10.util.ArrayList;
import x10x.vector.Point3d;
import au.edu.anu.chem.Atom;

/**
 * This class represents an atom for the purposes of
 * Quantum Chemistry simulations
 *
 * @author: V.Ganesh
 */
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
    
    var basisFunctions:ArrayList[ContractedGaussian];

    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian]) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() : ArrayList[ContractedGaussian] = basisFunctions;
}

