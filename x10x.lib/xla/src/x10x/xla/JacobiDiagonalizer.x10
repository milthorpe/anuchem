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
package x10x.xla;

import x10.compiler.Inline;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

/**
 * Based on JacobiDiagonalizer.java (org.meta.math.la.JacobiDiagonalizer)
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public class JacobiDiagonalizer {
    var eigenValuesVec:Vector; 
    var eigenVectorsMat:Matrix;

    val maxIterations : Int = 100;

    public def diagonalize(mat:Matrix) : void {
       val matrix = mat.getMatrix();
       val n  = mat.getRowCount();

       eigenVectorsMat = new Matrix(n);
       eigenVectorsMat.makeIdentity();

       val eigenVectors = eigenVectorsMat.getMatrix();

       // operate on a clone of the matrix
       val aMat = new Matrix(mat);
       val a = aMat.getMatrix();

       eigenValuesVec = new Vector(n); 
       val eigenValues = eigenValuesVec.getVector();
      
       val bVec = new Vector(n);
       val zVec = new Vector(n);
       val b = bVec.getVector();
       val z = zVec.getVector();

        for (i in 0..a.region.max(0)) {
            eigenValues(i) = b(i) = a(i,i);
            z(i) = 0.0;
        }

       var zeroTolerance:Double = 0.0;

       for(var sweepsIdx:Int = 0; sweepsIdx<maxIterations; sweepsIdx++) {
          val sum = aMat.sumOffDiagonal();
          if (sum == 0.0) break;  // if off diagonal elements are zero, we stop

          if (sweepsIdx < 3) zeroTolerance = 0.2 * sum / (n*n);
          else               zeroTolerance = 0.0;

          val sweeps = sweepsIdx;
          val zeroTol = zeroTolerance;

          for(var ip:Int = 0; ip<n-1; ip++) { for(var iq:Int = ip+1; iq<n; iq++) {
             if ((iq >= ip+1) && (ip < n-1)) { 

                val g:Double = 100.0 * Math.abs(a(ip,iq));

                if ((sweeps > 4) 
                    && (((Math.abs(eigenValues(ip))+g) as Double as Float) == ((Math.abs(eigenValues(ip))) as Double as Float))
                    && (((Math.abs(eigenValues(iq))+g) as Double as Float) == ((Math.abs(eigenValues(iq))) as Double as Float))) {
                   a(ip,iq) = 0.0;
                } else if (Math.abs(a(ip,iq)) > zeroTol) {
                   var h:Double = eigenValues(iq) - eigenValues(ip);
                   var t:Double, theta:Double;

                   if ((((Math.abs(h)+g)) as Double as Float) == (Math.abs(h) as Double as Float)) {
                     t = a(ip,iq) / h;
                   } else {
                     theta = 0.5 * h / a(ip,iq);
                     t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta*theta));
                     
                     if (theta < 0.0) t = -t;

                     val cos = 1.0 / Math.sqrt(1.0 + t*t);
                     val sin = t * cos;
                     val tau = sin / (1.0 + cos);
                     h   = t * a(ip, iq);

                     z(ip) -= h; z(iq) += h;
                     eigenValues(ip) -= h; eigenValues(iq) += h;

                     a(ip,iq) = 0.0;

                     for(var j:Int = 0; j<ip; j++)    doRotate(a, j, ip, j, iq, sin, tau);
                     for(var j:Int = ip+1; j<iq; j++) doRotate(a, ip, j, j, iq, sin, tau);
                     for(var j:Int = iq+1; j<n; j++)  doRotate(a, ip, j, iq, j, sin, tau);
                     for(var j:Int = 0; j<n; j++)     doRotate(eigenVectors, j, ip, j, iq, sin, tau);
                   } // end if
                } // end if
             } // end if
           }} // end for

          for([ip] in b) {
             b(ip) += z(ip);
             eigenValues(ip) = b(ip);
             z(ip) = 0.0;
          }
       } // end for

       sortEigenValues(); 
       eigenVectorsMat = eigenVectorsMat.transpose();
    }

    private @Inline static def doRotate(a:Array[Double](2){rect,zeroBased}, i:Int, j:Int, k:Int, l:Int,
                                sin:Double, tau:Double) : void {
       val g = a(i,j);
       val h = a(k,l);
       a(i,j) = g - sin * (h + g*tau);
       a(k,l) = h + sin * (g - h*tau);
    }

    private def sortEigenValues() {
        var i:Int, j:Int, k:Int;
        var p:Double;
        val n = eigenValuesVec.getSize();
        val eigenValues = eigenValuesVec.getVector();
        val eigenVectors = eigenVectorsMat.getMatrix();
        
        for(i=0; i<n; i++) {
            p = eigenValues(k=i);
            
            for(j=i+1; j<n; j++) 
                if (eigenValues(j) <= p) p = eigenValues(k=j);
            
            if (k!=i) { // swap
                eigenValues(k) = eigenValues(i);
                eigenValues(i) = p;
                
                for(j=0; j<n; j++) { // swap eigenVectors
                    p = eigenVectors(j,i);
                    eigenVectors(j,i) = eigenVectors(j,k);
                    eigenVectors(j,k) = p;
                } // end for
            } // end if
        } // end for
    }

    public def getEigenValues() : Vector = eigenValuesVec;

    public def getEigenVectors() : Matrix = eigenVectorsMat; 
}

