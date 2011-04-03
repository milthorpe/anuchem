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

import x10.util.ArrayList;

import x10x.vector.Point3d;

/**
 * This class represents a Ring structure in a Molecule
 *
 * @author: V.Ganesh
 */
public class Ring[T]{T <: Atom} {
    val atomList = new ArrayList[T](); 
    var planar:Boolean;

    public def this() { 
        planar = true;   
    }

    public def addAtom(atm:T) : void {
        atomList.add(atm); 
    }

    public def getAtom(index:Int) : T = atomList.get(index);
    public def getAtoms() = atomList;
    public def getSize() = atomList.size();

    public def isPlanar() = planar;
    public def setPlanar(p:Boolean) { planar = p; }

    public def toString() : String {
       var str:String = "[ ";

       for(atm:T in atomList)
         str += atm.toString() + " ";
       
       str += "]\n";

       return str;
    }
}

