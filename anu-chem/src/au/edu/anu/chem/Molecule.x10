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
import x10.compiler.Ifdef;
import x10.compiler.Ifndef;

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

    public def getCentreOfVDW(roZ:Double):Point3d {
        var step:Double = 0.1/roZ; // Atomic Unit
        val initialguess =  centreOfMass();
        val ai = AtomInfo.getInstance();
        var x:Double = -initialguess.i, y:Double = -initialguess.j, z:Double = -initialguess.k;
        var isConverged:Boolean=false;

        while (isConverged==false || step>.0001/roZ) {
            var rad1:Double=0.,rad2:Double=0.,dx:Double=0.,dy:Double=0.,dz:Double=0.;
            
            for(atm:T in atomList) {
                val atmvec= Vector3d(atm.centre.i+x,atm.centre.j+y,atm.centre.k+z);
                val distance=atmvec.magnitude()+ai.getVdwRadius(atm)/roZ;
                if (rad1<distance) { 
                    rad1=distance;
                    dx= -atm.centre.i/atmvec.magnitude()*step;
                    dy= -atm.centre.j/atmvec.magnitude()*step;
                    dz= -atm.centre.k/atmvec.magnitude()*step;
                }
            }
            for(atm:T in atomList) {
                val atmvec= Vector3d(atm.centre.i+x+dx,atm.centre.j+y+dy,atm.centre.k+z+dz);
                val distance=atmvec.magnitude()+ai.getVdwRadius(atm)/roZ;
                if (rad2<distance) rad2=distance;                    
            }
            @Ifdef("__DEBUG__") {Console.OUT.printf("rad1=%e, rad2=%e, step=%e\n",rad1,rad2,step);}
            if (rad2<rad1) {x+=dx; y+=dy;z+=dz; isConverged=false;}
            else {isConverged=true; step*=.1;}
        }
        return Point3d(-x,-y,-z);
    }

    /**
     * Translates and rotates this molecule to standard nuclear orientation.
     * @see Gill, P.M.W., Johnson, B.G. and Pople, J.A. (1993).
     *   "A standard grid for density functional calculations".
     *    Chem. Phys. Lett. 209 (5-6) pp.506-512  (Appendix A)
     */
    public def transformToSNO(roZ:Double) {
        // TODO use different method for molCentre according to input file specification.
        val molCentre = getCentreOfVDW(roZ); //getCentreOfNuclearCharge();
        val translation = Vector3d(-molCentre.i, -molCentre.j, -molCentre.k);
        if (translation.magnitude() > 1.0e-5/roZ) { // generally geometry files contain 3 decimal places
            Console.OUT.printf("\ntranslated molecule by [%.4g, %.4g, %.4g]\n\n", translation.i, translation.j, translation.k);
        }
        for(atm:T in atomList) {
            atm.centre += translation;
        }
        // TODO rotation
    }

    public def getRadius(roZ:Double):Double {
        var rad:Double=0.;
        val ai = AtomInfo.getInstance();
        for(atm:T in atomList) {
            val atmvec= Vector3d(atm.centre.i,atm.centre.j,atm.centre.k);
            val distance=atmvec.magnitude()+ai.getVdwRadius(atm)/roZ;
            if (rad<distance) rad=distance;
        }
        Console.OUT.printf("radius = %f  scaled radius = %f\n",rad*roZ,rad);
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

