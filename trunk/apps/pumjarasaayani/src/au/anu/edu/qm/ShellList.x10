/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
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
    global val shellList:HashMap[Int, Shell{self.at(this)}]{self.at(this)};
    global val shellPairs = new GrowableRail[Pair[Int,Int]]();
    
    var powerList:Rail[ValRail[Power]]{self.at(this)};

    var maxam:Int;

    public def this() { 
        shellList = new HashMap[Int, Shell{self.at(this)}]();
        maxam = 0;
    }

    public def initPowerList() : void {
        val maxam4 = (maxam*4)+2;
 
        val pl = new GrowableRail[ValRail[Power]](); 
        val pList = PowerList.getInstance() as PowerList{self.at(this)};
        for(var i:Int=0; i<=maxam4; i++)
           pl.add(pList.generatePowerList(i)); 

        powerList = pl.toRail();
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
    //public def getNumberOfShells() = shellList.size();

    public def getNumberOfShellPrimitives() : Int { 
         var n:Int = 0;
         for(shell in shellList.keySet()) 
             n += shellList.getOrElse(shell, null).getNumberOfShellPrimitives();
         return n;
    }

    public def getShellPrimitives() : GrowableRail[ContractedGaussian]! {
         val shellPrimitives = new GrowableRail[ContractedGaussian]();

         for(shell in shellList.keySet()) { 
            val sp = shellList.getOrElse(shell, null).getShellPrimitives();
            
            for(var i:Int=0; i<sp.length(); i++) shellPrimitives.add(sp(i));
         } // end for

         return shellPrimitives;
    } 

    var shellPairsGenerated:Boolean = false;

    public def getShellPairs() : GrowableRail[Pair[Int, Int]]! {
         if (!shellPairsGenerated) {
           val noOfShells = getNumberOfShellPrimitives();
           // Console.OUT.println("No of shells : " + noOfShells);
           for(var a:Int=0; a<noOfShells; a++) {
              for(var b:Int=0; b<noOfShells; b++) {
                 shellPairs.add(Pair[Int,Int](a,b));
              } // end for
           } // end for
           shellPairsGenerated = true;
         } // end if

         return shellPairs;
    }

    public def getShell(am:Int) = shellList.getOrElse(am, null);
    public def getPowers(am:Int) = powerList(am);
}

