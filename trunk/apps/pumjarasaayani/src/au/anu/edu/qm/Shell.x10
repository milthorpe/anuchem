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

import x10.util.*;

/**
 * Shell.x10
 *
 * Class representing a gaussian shell
 *
 * @author: V.Ganesh
 */
public class Shell { 
    val angularMomentum:Int;
    global val shellPrimitives:GrowableRail[ContractedGaussian]{self.at(this)};

    public def this(am:Int) { 
       angularMomentum = am;
       shellPrimitives = new GrowableRail[ContractedGaussian]();
    } 

    public def getAngularMomentum() = angularMomentum;
    public def getShellPrimitives() = shellPrimitives;
    public def getNumberOfShellPrimitives() = shellPrimitives.length();

    public def getShellPrimitive(i:Int) = shellPrimitives(i);

    public def addShellPrimitive(cg:ContractedGaussian) : void {
       shellPrimitives.add(cg);
    }
}

