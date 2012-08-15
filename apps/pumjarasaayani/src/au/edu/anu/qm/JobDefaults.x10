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

/**
 * SCF Job defaults
 *
 * @author V. Ganesh, milthorpe, T. Limpanuparb
 */
public class JobDefaults {
    public static SCF_METHOD_DIIS = "diis";
    public static SCF_METHOD_ROOTHAAN = "roothaan";

    public var scfMethod:String = SCF_METHOD_DIIS;

    public var maxIterations:Int;
    public var energyTolerance:Double;
    public var diisStartThreshold:Double;
    public var diisConvergenceThreshold:Double;
    public var diisSubspaceSize:Int;
    public var gMatrixParallelScheme:Int;
    public var useMta:Boolean;

    public static GUESS_CORE = "core";
    public static GUESS_SAD = "sad";
    public var guess:String;

// RO Var
    public var roOn:Int;
    public var compareRo:Boolean;
    public var roN:Int;
    public var roNK:Int;
    public var roL:Int;
    public var roZ:Double;

    public var omega:Double;
    public var thresh:Double;
    public var roThresh:Double;
    public var centering:Int; // See Q-CHEM Job description

    private def this() { 
        maxIterations = 100;
        energyTolerance = 1e-5;
        diisStartThreshold = 0.1;
        diisConvergenceThreshold = 1e-5; // Q-CHEM default
        diisSubspaceSize = 15; // Q-CHEM default
        gMatrixParallelScheme = GMatrix.DEFAULT_GMATTYPE;
        useMta = false;

        guess=JobDefaults.GUESS_SAD;
        
        roOn=0;
        compareRo = false;

        roN=10;
        roNK=-1; // see GmatrixRoMem.x10
        roL=10;
        roZ=1.0;

        omega = .1;
        thresh = 1e-8;
        roThresh = 1e-8;

    }

    private static _theInstance = new JobDefaults();
    public static def getInstance() = _theInstance;
}

