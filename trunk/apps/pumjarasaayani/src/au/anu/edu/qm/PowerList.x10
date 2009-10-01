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

    val powerList:HashMap[String,Array[Power]{rank==1}];
    var initCalled:Boolean;

    private def this() { 
       initCalled = false;
       powerList = new HashMap[String,Array[Power]{rank==1}](); 
    }

    private def init() : void {
       powerList.put("S", generatePowerList(0));
       powerList.put("P", generatePowerList(1));
       powerList.put("D", generatePowerList(2));
       powerList.put("F", generatePowerList(3)); 

       initCalled = true;
    }

    public def generatePowerList(maxAngularMomentum:Int) {
        var n:Int = 0;
        switch(maxAngularMomentum) {
          case 0:
               n = 1; break;
          case 1:
               n = 3; break;
          case 2:
               n = 6; break;
          case 3:
               n = 10; break;
        } // end switch .. case

        var pList:Array[Power]{rank==1};
        pList = Array.make[Power]([0..n]);
     
        var idx:Int = 0;
        for(var i:Int=maxAngularMomentum; i>=0; i--)
            for(var j:Int=maxAngularMomentum-i; j>=0; j--)
                pList(idx++) = new Power(i, j, maxAngularMomentum-i-j);

        return pList;
    }


    private static val _theInstance = new PowerList();

    public static def getInstance() { 
       if (!_theInstance.initCalled) _theInstance.init();

       return _theInstance;
    }

    public def getPowers(orbitalType:String) : Array[Power]{rank==1} {
       return (powerList.get(orbitalType) as Array[Power]{rank==1});
    }
}

