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
 *   <title>
 *   basis <6-31gdp>
 * + charge <n>
 * + multiplicity <n> [if charge specified]
 * + conversion <c> [optional unit conversion factor]
 *   molecule
 *   <symbol> <x> <y> <z>
 *   ...
 *   end
 * + scf_method [diis|roothaan]
 * + scf_max <n>
 * + scf_converge <x>
 * + diis_start <x>
 * + diis_converge <x>
 * + diis_subspace size <n>
 * + gmat_parallel_scheme <n>
 * + fragment mta
 * + scratch <filename>
 * + checkpoint <filename>
 * + guess [core|sad|file]
 * * print [hcore|overla|orthonorm|2e] [<filename>]
 * * scf_print [mo|density] [<filename>]
 * * scf_final_print [mo|density] [<filename>]
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
        // defaults
        val jd = JobDefaults.getInstance();
        var conversion:Double = 1.0;
        var charge:Int = 0;
        var multiplicity:Int = 1;
        jd.roOn=0;
        jd.roN=10;
        jd.roL=10;
        jd.roZ=1.0;
        jd.guess=0;

        if (line.startsWith("charge")) {
            charge = getIntParam(line);
            line = fil.readLine();
            val multWords = StringSplitter.splitOnWhitespace(line);
            if (!multWords(0).equals("multiplicity")) {
                throw new Exception("Invalid input: must specify multiplicity for charged molecule. Next line was:\n"+line);
            }
            multiplicity = Int.parseInt(multWords(1));
            line = fil.readLine();
        }

        if (line.startsWith("units")) {
            conversion = getDoubleParam(line);
            line = fil.readLine();
        }
 
        if (line.startsWith("RO_Z")) {
            jd.roZ = getDoubleParam(line); 
            line = fil.readLine();
        }
 
        if (!line.startsWith("molecule")) {
            throw new Exception("Invalid input: must specify molecule. Next line was:\n"+line);
        }

        molecule = new Molecule[QMAtom](title,charge,multiplicity);
        

        line = fil.readLine();
        try {
            while (line != null && !line.startsWith("end")) {
                val wrd = StringSplitter.splitOnWhitespace(line);

                molecule.addAtom(new QMAtom(wrd(0), 
                                   Point3d(Double.parseDouble(wrd(1))/conversion/jd.roZ,
                                            Double.parseDouble(wrd(2))/conversion/jd.roZ,
                                            Double.parseDouble(wrd(3))/conversion/jd.roZ
                                    )
                                ));

                line = fil.readLine();
            }
            line = fil.readLine();
            while (line != null && !line.startsWith("end")) {
                if (line.startsWith("scf_max")) {
                    jd.maxIterations = getIntParam(line);
                } else if (line.startsWith("scf_converge")) {
                    jd.energyTolerance = getDoubleParam(line);
                } else if (line.startsWith("diis_start")) {
                    jd.diisStartThreshold = getDoubleParam(line);
                } else if (line.startsWith("diis_converge")) {
                    jd.diisConvergenceThreshold = getDoubleParam(line);
                } else if (line.startsWith("diis_subspace_size")) {
                    jd.diisSubspaceSize = getIntParam(line);
                } else if (line.startsWith("gmat_parallel_scheme")) {
                    jd.gMatrixParallelScheme = getIntParam(line);
                } else if (line.startsWith("fragment mta")) {
                    jd.useMta = true;
                } else if (line.startsWith("RO_N")) {
                    jd.roN = getIntParam(line);
                } else if (line.startsWith("RO_L")) {
                    jd.roL = getIntParam(line);
                } else if (line.startsWith("Center")) {
                    jd.centering = getIntParam(line);
                } else if (line.startsWith("USE_RO")) {
                    jd.roOn = getIntParam(line);
                } else if (line.startsWith("GUESS")) {
                    jd.guess = getIntParam(line);
                }

                line = fil.readLine();
            }
        } catch (e:EOFException) {
            // no more atoms / job directives
        }

        fil.close();
    }

    private static def getIntParam(line:String) = Int.parseInt(StringSplitter.splitOnWhitespace(line)(1));
    private static def getDoubleParam(line:String) = Double.parseDouble(StringSplitter.splitOnWhitespace(line)(1));

    public def getMolecule():Molecule[QMAtom] = molecule;
    public def getBasisName():String = basisName;
}

