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
package au.edu.anu.qm;

import x10.compiler.Inline;
import x10.util.ArrayList;
import x10.util.HashMap;
import au.edu.anu.chem.Molecule;

/**
 * Structure storing Shell list
 *
 * @author: V.Ganesh, milthorpe
 */
public struct ShellList { 
    private val shellList:HashMap[Int, Shell];
    private val powerList:Rail[Rail[Power]];
    private val maxam:Int;
    private val maxDcon:Int; // Max degree of contraction

    public def this(molecule:Molecule[QMAtom]) {
        // init shell list
        val shellList = new HashMap[Int, Shell]();
        var maxam:Int = 0n; var maxDcon:Int = 0n;
        for(var atmno:Int=0n; atmno<molecule.getNumberOfAtoms(); atmno++) {
            val atom = molecule.getAtom(atmno);
            val bfs  = atom.getBasisFunctions();
            val nbf  = bfs.size();

            for(var i:Int=0n; i<nbf; i++) {
                val cg = bfs.get(i);
                val am = cg.getMaximumAngularMomentum();
                val dcon = cg.getMaximumDegreeOfContraction();
                maxam  = Math.max(am, maxam);
                maxDcon = Math.max(dcon, maxDcon);

                var shell:Shell = shellList.getOrElse(am, null);
                if (shell == null) {
                   shell = new Shell(am);
                   shellList.put(am, shell);
                }
                shell.addShellPrimitive(cg);
            } // end for
        } // end for

        // init power list
        val maxam4 = (maxam*4n)+2n;
        val powerList = new Rail[Rail[Power]](maxam4+1n); 
        val pList = PowerList.getInstance();
        for(i in 0n..maxam4) powerList(i) = pList.generatePowerList(i); 

        this.maxam = maxam;
        this.maxDcon = maxDcon;
        this.shellList = shellList;
        this.powerList = powerList;
    }

    public def getNumberOfShellPrimitives() : Int { 
        var n:Int = 0n;
        for(shell in shellList.keySet()) 
            n += shellList.getOrElse(shell, null).getNumberOfShellPrimitives();
        return n;
    }

    public def getShellPrimitives() : Rail[ContractedGaussian] {
        val shellPrimitives = new ArrayList[ContractedGaussian]();

        for(shell in shellList.keySet()) { 
            val sp = shellList.getOrElse(shell, null).getShellPrimitives();
            for(var i:Int=0n; i<sp.size(); i++) shellPrimitives.add(sp(i));
        }

        return shellPrimitives.toRail();
    }

    public def getNumberOfShellPairs() : Int {
        val nPrimitives = getNumberOfShellPrimitives();
        return nPrimitives * nPrimitives;
    }

    public @Inline def getMaximumAngularMomentum() = maxam;
    public @Inline def getMaximumDegreeOfContraction() = maxDcon;
    public @Inline def getShell(am:Int) = shellList.getOrElse(am, null);
    public @Inline def getPowers(am:Int) = powerList(am);
}

