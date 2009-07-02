package au.anu.edu.qm;

import x10.util.*;

public class PowerList { 

    var powerList:HashMap[String,Array[Power]{rank==1}];
    var initCalled:Boolean;

    private def this() { 
       initCalled = false;
    }

    private def init() : void {
       powerList = new HashMap[String,Array[Power]{rank==1}](); 

       val sList:Array[Power]{rank==1} = sList = Array.make[Power]([0..1]);
       sList(0) = new Power(0,0,0);
       powerList.put("S", sList);

       val pList:Array[Power]{rank==1} = Array.make[Power]([0..3]);
       pList(0) = new Power(1,0,0);
       pList(1) = new Power(0,1,0);
       pList(2) = new Power(0,0,1);
       powerList.put("P", pList);

       val dList:Array[Power]{rank==1} = Array.make[Power]([0..6]);
       dList(0) = new Power(2,0,0);
       dList(1) = new Power(0,2,0);
       dList(2) = new Power(0,0,2);
       dList(3) = new Power(1,1,0);
       dList(4) = new Power(0,1,1);
       dList(5) = new Power(1,0,1);
       powerList.put("D", dList);

       val fList:Array[Power]{rank==1} = Array.make[Power]([0..10]);
       fList(0) = new Power(3,0,0);
       fList(1) = new Power(2,1,0);
       fList(2) = new Power(2,0,1);
       fList(3) = new Power(1,2,0);
       fList(4) = new Power(1,1,1);
       fList(5) = new Power(1,0,2);
       fList(6) = new Power(0,3,0);
       fList(7) = new Power(0,2,1);
       fList(8) = new Power(0,1,2);
       fList(9) = new Power(0,0,3);
       powerList.put("F", fList); 

       initCalled = true;
    }

    static val _theInstance:PowerList = new PowerList();

    public static def getInstance() : PowerList { 
       if (!_theInstance.initCalled) _theInstance.init();

       return _theInstance;
    }

    public def getPowers(orbitalType:String) : Array[Power]{rank==1} {
       return (powerList.get(orbitalType) as Array[Power]{rank==1});
    }
}

