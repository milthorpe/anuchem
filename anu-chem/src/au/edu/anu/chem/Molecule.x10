/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.chem;

import x10.util.Pair;
import x10.util.ArrayList;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents a Molecule
 * for the purposes of computational chemistry codes.
 *
 * @author: V.Ganesh
 */
public class Molecule[T]{T <: Atom} {
    val atomList = new ArrayList[T](); 
    val name : String;
    val charge : Int;
    val multiplicity : Int;

    val ringList = new ArrayList[Ring[T]]();

    /** 
     * Measures the maximum absolute value of any coordinate x,y,z
     * of all atoms. This is used to estimate a rough cubic box size.
     */
    private var maxExtent : Double = 0.0;
    // Require for Fragment.x10:35-38:
    public def this() { 
        name = "unknown"; 
        this.charge = 0; 
        this.multiplicity = 1;
    }

    public def this(name:String) {
        this(name, 0, 1);
    }

    public def this(name:String, c:Int, m:Int) {
        this.name = name;
        this.charge = c;
        this.multiplicity = m;
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
    public def getCharge() : Int = charge;
    public def getMultiplicity() : Int = multiplicity;

    public def addRing(r:Ring[T]) { ringList.add(r); }
    public def getRings() : ArrayList[Ring[T]] = ringList;

    public def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:T in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne-charge;
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

    public def getCentreOfNuclearCharge():Point3d {
        val ai = AtomInfo.getInstance();
        var x:Double = 0.0, y:Double = 0.0, z:Double = 0.0;
        var chargeSum:Double = 0.0;

        for(atm:T in atomList) {
            val nuclearCharge = ai.getAtomicNumber(atm);
            x += nuclearCharge * atm.centre.i;
            y += nuclearCharge * atm.centre.j;
            z += nuclearCharge * atm.centre.k;

            chargeSum += nuclearCharge;
        }

        return Point3d(x/chargeSum,y/chargeSum,z/chargeSum);
    }

    /**
     * Translates and rotates this molecule to standard nuclear orientation.
     * @see Gill, P.M.W., Johnson, B.G. and Pople, J.A. (1993).
     *   "A standard grid for density functional calculations".
     *    Chem. Phys. Lett. 209 (5-6) pp.506-512  (Appendix A)
     */
    public def transformToSNO() {
        val zeroMoment = getCentreOfNuclearCharge();
        val translation = Vector3d(-zeroMoment.i, -zeroMoment.j, -zeroMoment.k);
        if (translation.magnitude() > 1.0e-12) {
            Console.OUT.printf("\ntranslated molecule by [%.4g, %.4g, %.4g]\n\n", translation.i, translation.j, translation.k);
        }
        for(atm:T in atomList) {
            atm.centre += translation;
        }
        // TODO rotation
    }

    public def getRadius():Double {
        var rad:Double=0.;
        val ai = AtomInfo.getInstance();
        for(atm:T in atomList) {
            val distance=atm.centre.magnitude()+ai.getVdwRadius(atm.symbol);
            if (rad<distance) rad=distance;
        }
        return rad; 
    }

    public def toString() : String {
        var str:String = "";

        if (charge != 0) {
            str += "charge " + charge + "\n"
                + "multiplicity " + multiplicity + "\n";
        }

        for(atm:T in atomList)
         str += atm.toString() + "\n";

        return str;
    }
}

