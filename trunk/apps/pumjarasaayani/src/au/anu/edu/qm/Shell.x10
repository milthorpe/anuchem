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
    global val shellPrimitives:ArrayList[ContractedGaussian]{self.at(this)};

    public def this(am:Int) { 
       angularMomentum = am;
       shellPrimitives = new ArrayList[ContractedGaussian]();
    } 

    public def getAngularMomentum() = angularMomentum;
    public def getShellPrimitives() = shellPrimitives;
    public def getNumberOfShellPrimitives() = shellPrimitives.size();

    public def addShellPrimitive(cg:ContractedGaussian) : void {
       shellPrimitives.add(cg);
    }
}

