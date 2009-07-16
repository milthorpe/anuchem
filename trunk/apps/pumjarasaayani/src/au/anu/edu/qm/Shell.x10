package au.anu.edu.qm;

import x10.util.*;

public class Shell { 
    var angularMomentum:Int;
    var shellPrimitives:ArrayList[ContractedGaussian];

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

