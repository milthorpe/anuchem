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
package au.edu.anu.chem;

import x10.util.Pair;
import x10.util.ArrayList;

import x10x.vector.Point3d;

/**
 * This class represents a Molecule
 * for the purposes of computational chemistry codes.
 *
 * @author: V.Ganesh
 */
public class Molecule[T]{T <: Atom} {
    val atomList = new ArrayList[T](); 
    val name : String;

    val ringList = new ArrayList[Ring[T]]();

    /** 
     * Measures the maximum absolute value of any coordinate x,y,z
     * of all atoms. This is used to estimate a rough cubic box size.
     */
    private var maxExtent : Double = 0.0;

    public def this() { 
        name = "unknown";
    }

    public def this(name:String) {
        this.name = name;
    }

    public def getName() = this.name;

    public def addAtom(atm:T) : void {
        atomList.add(atm); 
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.i));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.j));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.k));      
    }

    public def getAtom(index:Int) : T = atomList.get(index);
    public def getAtoms() = atomList;
    public def getNumberOfAtoms() : Int = atomList.size();

    public def addRing(r:Ring[T]) { ringList.add(r); }
    public def getRings() : ArrayList[Ring[T]] = ringList;

    public def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:T in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public def getMaxExtent() = maxExtent;

    public def getCoords() : Rail[Pair[String,Point3d]] {
        val coords = new ArrayList[Pair[String,Point3d]](atomList.size());
        for(atom in atomList) {
            coords.add(Pair[String,Point3d](atom.symbol, atom.centre));
        }
        return coords.toArray();
    }

    public def centreOfMass() : Point3d {
        val ai = AtomInfo.getInstance();
        var x:Double = 0.0, y:Double = 0.0, z:Double = 0.0;
        var massSum:Double = 0.0;

        for(atm:T in atomList) {
            val mass = ai.getAtomicMass(atm);
            x += mass * atm.centre.i;
            y += mass * atm.centre.j;
            z += mass * atm.centre.k;

            massSum += mass;
        }

        return Point3d(x/massSum,y/massSum,z/massSum);
    }

    public def toString() : String {
       var str:String = "";

       for(atm:T in atomList)
         str += atm.toString() + "\n";

       return str;
    }
}

