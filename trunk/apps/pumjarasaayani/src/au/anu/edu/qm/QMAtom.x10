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

import x10.util.ArrayList;
import x10x.vector.Point3d;
import au.edu.anu.chem.Atom;

/**
 * QMAtom.x10
 *
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
    
    var basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;
}

