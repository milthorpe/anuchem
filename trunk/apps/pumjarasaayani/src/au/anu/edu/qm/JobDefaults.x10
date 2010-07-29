/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

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

