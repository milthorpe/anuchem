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

import x10.io.File;
import x10.io.FileReader;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.StringSplitter;
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
       Console.OUT.println("NOTE: the input file is expected to be in ** a.u. ** units");

       val fil = new FileReader(new File(inpFile));

       val noOfAtoms = Int.parseInt(fil.readLine());

       val words = StringSplitter.splitOnWhitespace(fil.readLine());
       basisName = words(1);
       val charge:Int;
       val multiplicity:Int;
       val conversion:Double;
       if (words.size > 2) charge = Int.parseInt(words(2)); else charge = 0;
       if (words.size > 3) multiplicity = Int.parseInt(words(3)); else multiplicity = 1;
       if (words.size > 4) conversion = Double.parseDouble(words(4)); else conversion=1.0;
       molecule = new Molecule[QMAtom](words(0),charge,multiplicity);    

       for(var i:Int=0; i<noOfAtoms; i++) { 
         val wrd = StringSplitter.splitOnWhitespace(fil.readLine());

         molecule.addAtom(new QMAtom(wrd(0), 
                                       Point3d(Double.parseDouble(wrd(1))/conversion,
                                               Double.parseDouble(wrd(2))/conversion,
                                               Double.parseDouble(wrd(3))/conversion
                                    )
                            ));
       }

       fil.close();
    }

    public def getMolecule()  : Molecule[QMAtom] = molecule;
    public def getBasisName() : String   = basisName;
}

