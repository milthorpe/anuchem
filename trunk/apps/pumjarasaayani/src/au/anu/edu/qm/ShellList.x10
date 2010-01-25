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
    global var shellList:HashMap[Int, Shell{self.at(this)}]{self.at(this)};

    public def this() { 
    }

    public def make() { 
        shellList = new HashMap[Int, Shell{self.at(this)}]();
    } 

    public def addShellPrimitive(cg:ContractedGaussian{self.at(this)}) : void {
        var shell:Shell{self.at(this)} = getShell(cg.getTotalAngularMomentum());
        if (shell == null) {
           shell = new Shell(cg.getTotalAngularMomentum());
           shellList.put(cg.getTotalAngularMomentum(), shell);
        }
        shell.addShellPrimitive(cg);
    }

    public def getNumberOfShells() : Int = shellList.size();
    public def getShell(am:Int) = shellList.getOrElse(am, null);
}

