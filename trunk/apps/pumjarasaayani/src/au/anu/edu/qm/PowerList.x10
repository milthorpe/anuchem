/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

import x10.util.*;

/**
 * PowerList.x10
 *
 * Represents gaussian powers terms for various orbital types
 *
 * @author: V.Ganesh
 */
public class PowerList { 

    val powerList = new HashMap[String,ValRail[Power]]() as HashMap[String,ValRail[Power]]!; 

    private def this() {
       powerList.put("S", generatePowerList(0));
       powerList.put("P", generatePowerList(1));
       powerList.put("D", generatePowerList(2));
       powerList.put("F", generatePowerList(3)); 
    }

    public def generatePowerList(maxAngularMomentum:Int) : ValRail[Power] {
        var n:Int = ((maxAngularMomentum+1)*(maxAngularMomentum+2)/2);

        val pList = Rail.make[Power](n);

        var idx:Int = 0;
        // for(var i:Int=maxAngularMomentum; i>=0; i--) {
        //    for(var j:Int=maxAngularMomentum-i; j>=0; j--) {
        for(var i:Int=0; i<=maxAngularMomentum; i++) {
            for(var j:Int=0; j<=maxAngularMomentum-i; j++) {
                pList(idx++) = Power(i, j, maxAngularMomentum-i-j);
            }
        }

        return pList as ValRail[Power];
    }

    private static val _theInstance = new PowerList() as PowerList!;

    public static def getInstance() : PowerList! { 
       // return _theInstance as PowerList!;
       return new PowerList();
    }

    public def getPowers(orbitalType:String) = powerList.getOrElse(orbitalType, null);
}

