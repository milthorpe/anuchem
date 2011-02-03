/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.anu.edu.qm;

import x10.util.*;

/**
 * ShellList.x10
 *
 * Structure storing Shell list
 *
 * @author: V.Ganesh
 */
public class ShellList { 
    val shellList:HashMap[Int, Shell];
    
    var powerList:Array[Array[Power](1){rect,rail}](1){rect,rail};

    var maxam:Int;

    public def this() { 
        shellList = new HashMap[Int, Shell]();
        maxam = 0;
    }

    public def initPowerList() : void {
        val maxam4 = (maxam*4)+2;
 
        powerList = new Array[Array[Power](1){rect,rail}](maxam4); 
        val pList = PowerList.getInstance();
        for(var i:Int=0; i<=maxam4; i++)
           powerList(i) = pList.generatePowerList(i); 

    }

    public def addShellPrimitive(cg:ContractedGaussian) : void {
        val am = cg.getMaximumAngularMomentum();
        maxam  = Math.max(am, maxam);

        var shell:Shell = getShell(am);
        if (shell == null) {
           shell = new Shell(am);
           shellList.put(am, shell);
        }
        shell.addShellPrimitive(cg);
    }

    public def getMaximumAngularMomentum() = maxam;
    //public def getNumberOfShells() = shellList.size();

    public def getNumberOfShellPrimitives() : Int { 
         var n:Int = 0;
         for(shell in shellList.keySet()) 
             n += shellList.getOrElse(shell, null).getNumberOfShellPrimitives();
         return n;
    }

    public def getShellPrimitives() : ArrayList[ContractedGaussian] {
         val shellPrimitives = new ArrayList[ContractedGaussian]();

         for(shell in shellList.keySet()) { 
            val sp = shellList.getOrElse(shell, null).getShellPrimitives();
            
            for(var i:Int=0; i<sp.size(); i++) shellPrimitives.add(sp(i));
         } // end for

         return shellPrimitives;
    }

    public def getNumberOfShellPairs() : Int {
        val nPrimitives = getNumberOfShellPrimitives();
        return nPrimitives * nPrimitives;
    }

    public def getShell(am:Int) = shellList.getOrElse(am, null);
    public def getPowers(am:Int) = powerList(am);
}

