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
package au.edu.anu.qm;

import x10.compiler.NonEscaping;
import x10.util.HashMap;

/**
 * Represents gaussian powers terms for various orbital types
 *
 * @author: V.Ganesh
 */
public class PowerList { 

    val powerList = new HashMap[String,Rail[Power]](); 

    private def this() {
       powerList.put("S", generatePowerList(0n));
       powerList.put("P", generatePowerList(1n));
       powerList.put("D", generatePowerList(2n));
       powerList.put("F", generatePowerList(3n)); 
       powerList.put("G", generatePowerList(4n)); 
       powerList.put("H", generatePowerList(5n)); 
    }

    @NonEscaping 
    public final def generatePowerList(maxAngularMomentum:Int) : Rail[Power] {
        var n:Int = ((maxAngularMomentum+1n)*(maxAngularMomentum+2n)/2n);

        val pList = new Rail[Power](n);

        var idx:Long = 0;
        // for(var i:Int=maxAngularMomentum; i>=0; i--) {
        //    for(var j:Int=maxAngularMomentum-i; j>=0; j--) {
        for(var i:Int=0n; i<=maxAngularMomentum; i++) {
            for(var j:Int=0n; j<=maxAngularMomentum-i; j++) {
                pList(idx++) = Power(i, j, maxAngularMomentum-i-j);
            }
        }

        return pList;
    }

    private static val _theInstance = new PowerList();

    public static def getInstance() : PowerList { 
       // return _theInstance as PowerList!;
       return new PowerList();
    }

    public def getPowers(orbitalType:String) = powerList.getOrElse(orbitalType, null);
}

