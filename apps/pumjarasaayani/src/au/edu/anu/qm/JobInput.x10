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
import x10.io.EOFException;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.StringSplitter;
import x10x.vector.Point3d;

/**
 * This class parses an input file in the format described below.
 * Coordinates are assumed to be in a.u. if no conversion factor is supplied.
 *
 *  <title>
 *  basis 6-31gdp
 *  charge 0 [optional]
 *  multiplicity 1 [if charge specified]
 *  conversion 1.0 [optional unit conversion factor]
 *  molecule
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
        val fil = new FileReader(new File(inpFile));

        val title = fil.readLine();

        var line:String = fil.readLine();
        val basisWords = StringSplitter.splitOnWhitespace(line);
        if (!basisWords(0).equals("basis")) {
            throw new Exception("Invalid input: must specify basis. Next line was:\n"+line);
        }
        basisName = basisWords(1);
        line = fil.readLine();

        var charge:Int = 0;
        var multiplicity:Int = 1;
        if (line.startsWith("charge")) {
            charge = Int.parseInt(StringSplitter.splitOnWhitespace(line)(1));
            line = fil.readLine();
            val multWords = StringSplitter.splitOnWhitespace(line);
            if (!multWords(0).equals("multiplicity")) {
                throw new Exception("Invalid input: must specify multiplicity for charged molecule. Next line was:\n"+line);
            }
            multiplicity = Int.parseInt(multWords(1));
            line = fil.readLine();
        }

        var conversion:Double = 1.0;
        if (line.startsWith("units")) {
            conversion = Double.parseDouble(StringSplitter.splitOnWhitespace(line)(1));
            line = fil.readLine();
        }

        if (!line.startsWith("molecule")) {
            throw new Exception("Invalid input: must specify molecule. Next line was:\n"+line);
        }

        molecule = new Molecule[QMAtom](title,charge,multiplicity);    

        line = fil.readLine();
        while (line != null && line.trim().length() > 0) {
            val wrd = StringSplitter.splitOnWhitespace(line);

            molecule.addAtom(new QMAtom(wrd(0), 
                               Point3d(Double.parseDouble(wrd(1))/conversion,
                                        Double.parseDouble(wrd(2))/conversion,
                                        Double.parseDouble(wrd(3))/conversion
                                )
                            ));
            try {
                line = fil.readLine();
            } catch (e:EOFException) {
                // no more atoms
                break;
            }
        }

        fil.close();
    }

    public def getMolecule():Molecule[QMAtom] = molecule;
    public def getBasisName():String = basisName;
}

