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
package au.edu.anu.qm;

/**
 * SCF Job defaults
 *
 * @author V. Ganesh
 */
public class JobDefaults {
   private var maxIterations:Int;
   private var energyTolerance:Double;
   private var densityTolerance:Double; 

   private def this() { 
       maxIterations    = 100;
       energyTolerance  = 1e-5;
       densityTolerance = 1e-4;
   }

   private static _theInstance = new JobDefaults();
   public static def getInstance() = _theInstance;

   public def getMaxIterations() = maxIterations;
   public def setMaxIterations(it:Int) {
      maxIterations = maxIterations;
   }

   public def getEnergyTolerance() = energyTolerance;
   public def setEnergyTolerance(et:Double) {
      energyTolerance = et;
   }

   public def getDensityTolerance() = energyTolerance;
   public def setDensityTolerance(dt:Double) {
      densityTolerance = dt;
   }
}

