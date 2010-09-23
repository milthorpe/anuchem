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

import x10.io.*;
import au.edu.anu.chem.Molecule;
import x10x.vector.Point3d;

/**
 * This class expects the input file to be in following format.
 * Coordinates are expected to be in a.u.
 *
 *  <number of atoms>
 *  <title> <basis>
 *  symbol x y z
 *  ...
 */
public class JobInput { 
    var molecule:Molecule[QMAtom];
    var basisName:String;

    public def make(inpFile:String) { 
       readInp(inpFile);
    } 

    private def readInp(inpFile:String) { 
       x10.io.Console.OUT.println("NOTE: the input file is expected to be in ** a.u. ** units");

       val fil = new FileReader(new File(inpFile));

       val noOfAtoms = Int.parseInt(fil.readLine());
       val words = fil.readLine().split(" ");
 
       molecule = new Molecule[QMAtom](words(0));
       basisName = words(1);

       for(var i:Int=0; i<noOfAtoms; i++) { 
         val wrd = fil.readLine().split(" ");         

         molecule.addAtom(new QMAtom(wrd(0), 
                                       Point3d(Double.parseDouble(wrd(1)),
                                               Double.parseDouble(wrd(2)),
                                               Double.parseDouble(wrd(3))
                                    )
                            ));
       } // end while              

       fil.close();
    }

    public def getMolecule()  : Molecule[QMAtom] = molecule;
    public def getBasisName() : String   = basisName;
}

