/**
 * PowerList.x10
 *
 * Represents gaussian powers terms for various orbital types
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class PowerList { 

    global var powerList:HashMap[String,Array[Power{self.at(this)}]{self.at(this),rank==1}]{self.at(this)};
    global var initCalled:Boolean;

    private def this() {
       initCalled = false;
    }

    private def init() : void {
       powerList = new HashMap[String,Array[Power{self.at(this)}]{self.at(this),rank==1}](); 
       powerList.put("S", generatePowerList(0));
       powerList.put("P", generatePowerList(1));
       powerList.put("D", generatePowerList(2));
       powerList.put("F", generatePowerList(3)); 

       initCalled = true;
    }

    public def generatePowerList(maxAngularMomentum:Int) : Array[Power{self.at(this)}]{self.at(this),rank==1} {
        var n:Int = ((maxAngularMomentum+1)*(maxAngularMomentum+2)/2);

        val pList:Array[Power{self.at(this)}]{self.at(this),rank==1} = Array.make[Power{self.at(this)}]([0..n]);

        var idx:Int = 0;
        // for(var i:Int=maxAngularMomentum; i>=0; i--) {
        //    for(var j:Int=maxAngularMomentum-i; j>=0; j--) {
        for(var i:Int=0; i<=maxAngularMomentum; i++) {
            for(var j:Int=0; j<=maxAngularMomentum-i; j++) {
                var pwr:Power{self.at(this)} = new Power(i, j, maxAngularMomentum-i-j);
                pList(idx++) = pwr;
            }
        }

        return pList;
    }


    private global static val _theInstance = new PowerList();

    public static def getInstance() : PowerList { 
       if (!_theInstance.initCalled) _theInstance.init();

       return _theInstance;
    }

    public def getPowers(orbitalType:String) = powerList.getOrElse(orbitalType, null);
}

