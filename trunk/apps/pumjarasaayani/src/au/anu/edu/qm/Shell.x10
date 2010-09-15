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
    val shellPrimitives:GrowableRail[ContractedGaussian];

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

