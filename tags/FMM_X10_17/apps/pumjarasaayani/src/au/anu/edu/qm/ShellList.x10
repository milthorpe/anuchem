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
    val shellList:HashMap[Int, Shell];

    public def this() { 
        shellList = new HashMap[Int, Shell]();
    } 

    public def addShellPrimitive(cg:ContractedGaussian) : void {
        var shell:Shell = null; 
        
        try {
           shell = shellList.get(cg.getTotalAngularMomentum()) as Shell;
        } catch (e:NullPointerException) {
           shell = new Shell(cg.getTotalAngularMomentum());
           shellList.put(cg.getTotalAngularMomentum(), shell);
        } // end try .. catch 

        shell.addShellPrimitive(cg);
    }

    public def getNumberOfShell() : Int = shellList.size();
    public def getShell(am:Int) : Shell = shellList.get(am) as Shell;
}
