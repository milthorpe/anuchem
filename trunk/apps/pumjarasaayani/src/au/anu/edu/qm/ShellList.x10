/**
 * ShellList.x10
 *
 * Structure storing Shell list
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class ShellList { 
    global val shellList:HashMap[Int, Shell{self.at(this)}]{self.at(this)};
    global val powerList:ArrayList[ValRail[Power]]{self.at(this)};

    var maxam:Int;

    public def this() { 
        shellList = new HashMap[Int, Shell{self.at(this)}]();
        powerList = new ArrayList[ValRail[Power]]();
        maxam = 0;
    }

    public def initPowerList() : void {
        val maxam4 = (maxam*4)+2;
  
        val pList = PowerList.getInstance() as PowerList{self.at(this)};
        for(var i:Int=0; i<=maxam4; i++)
           powerList.add(pList.generatePowerList(i)); 
    }

    public def addShellPrimitive(cg:ContractedGaussian{self.at(this)}) : void {
        val am = cg.getMaximumAngularMomentum();
        maxam  = Math.max(am, maxam);

        var shell:Shell{self.at(this)} = getShell(am);
        if (shell == null) {
           shell = new Shell(am);
           shellList.put(am, shell);
        }
        shell.addShellPrimitive(cg);
    }

    public def getMaximumAngularMomentum() = maxam;
    public def getNumberOfShells() = shellList.size();
    public def getShell(am:Int) = shellList.getOrElse(am, null);
    public def getPowers(am:Int) = powerList.get(am);
}

