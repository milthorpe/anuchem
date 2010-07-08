/**
 * Shell.x10
 *
 * Class representing a gaussian shell 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

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

