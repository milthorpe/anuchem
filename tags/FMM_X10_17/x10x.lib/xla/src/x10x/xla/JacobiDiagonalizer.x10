package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

/**
 * Based on JacobiDiagonalizer.java (org.meta.math.la.JacobiDiagonalizer)
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public class JacobiDiagonalizer implements Diagonalizer {
    var eigenValuesVec : Vector; 
    var eigenVectorsMat : Matrix;
    var eigenValues  : Array[Double]{rank==1};
    var eigenVectors : Array[Double]{rank==2};
    var cos:Double, sin:Double, tau:Double;

    val maxIterations : Int = 100;

    public def diagonalize(mat:Matrix) : void {
       val matrix = mat.getMatrix();
       val n:Int = mat.getRowCount();

       eigenVectorsMat = Matrix.make(mat.dist());
       eigenVectorsMat.makeIdentity();

       eigenVectors = eigenVectorsMat.getMatrix();

       // clone the matrix, do not tamper the actual matrix
       val aMat = Matrix.make(mat.dist());
       val a = aMat.getMatrix();;
       ateach(var(i,j) in a.dist) a(i, j) = matrix(i, j);

       eigenValuesVec = Vector.make(n); 
       eigenValues = eigenValuesVec.getVector();
      
       val bVec = Vector.make(n);
       val zVec = Vector.make(n);
       val b = bVec.getVector();
       val z = zVec.getVector();

       for(var i:Int = 0; i<n; i++) { 
          eigenValues(i) = b(i) = a(i,i);
          z(i) = 0.0;
       } // end for

       var zeroTolerance:Double = 0.0;

       for(var sweeps:Int = 0; sweeps<maxIterations; sweeps++) {
          // sum off diagonal elements of A
          val sum = aMat.sumOffDiagonal();
          
          if (sum == 0.0) break;  // if off diagonal elements are zero, we stop

          if (sweeps < 3) zeroTolerance = 0.2 * sum / (n*n);
          else            zeroTolerance = 0.0;

          for(var ip:Int = 0; ip<n-1; ip++) {
             for(var iq:Int =ip+1; iq<n; iq++) {
                val g:Double = 100.0 * Math.abs(a(ip,iq));

                if ((sweeps > 4) 
                    && (((Math.abs(eigenValues(ip))+g) as Double as Float) == ((Math.abs(eigenValues(ip))) as Double as Float))
                    && (((Math.abs(eigenValues(iq))+g) as Double as Float) == ((Math.abs(eigenValues(iq))) as Double as Float))) {
                   a(ip,iq) = 0.0;
                } else if (Math.abs(a(ip,iq)) > zeroTolerance) {
                   var h:Double = eigenValues(iq) - eigenValues(ip);
                   var t:Double, theta:Double;

                   if ((((Math.abs(h)+g)) as Double as Float) == (Math.abs(h) as Double as Float)) {
                     t = a(ip,iq) / h;
                   } else {
                     theta = 0.5 * h / a(ip,iq);
                     t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta*theta));
                     
                     if (theta < 0.0) t = -t;

                     cos = 1.0 / Math.sqrt(1.0 + t*t);
                     sin = t * cos;
                     tau = sin / (1.0 + cos);
                     h   = t * a(ip, iq);

                     z(ip) -= h; z(iq) += h;
                     eigenValues(ip) -= h; eigenValues(iq) += h;

                     a(ip,iq) = 0.0;

                     for(var j:Int = 0; j<ip; j++)    doRotate(a, j, ip, j, iq);
                     for(var j:Int = ip+1; j<iq; j++) doRotate(a, ip, j, j, iq);
                     for(var j:Int = iq+1; j<n; j++)  doRotate(a, ip, j, iq, j);
                     for(var j:Int = 0; j<n; j++)     doRotate(eigenVectors, j, ip, j, iq);
                   } // end if
                } // end if
             } // end for
          } // end for          

          for(var ip:Int = 0; ip<n; ip++) {
             b(ip) += z(ip);
             eigenValues(ip) = b(ip);
             z(ip) = 0.0;
          } // end for
       } // end for

       sortEigenValues(); 

       var temp:Double;
       for(var i:Int = 0; i<n; i++) 
          for(var j:Int = i+1; j<n; j++) {
             temp = eigenVectors(i, j);
             eigenVectors(i, j) = eigenVectors(j, i);
             eigenVectors(j, i) = temp;
          }      
    }

    private def sortEigenValues() {
        var i:Int, j:Int, k:Int;
        var p:Double;
        val n = eigenValuesVec.getSize();
        
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

    private def doRotate(a:Array[Double]{rank==2}, i:Int, j:Int, k:Int, l:Int) : void {
       var g:Double = a(i,j);
       var h:Double = a(k,l);
 
       a(i,j) = g - sin * (h + g*tau);
       a(k,l) = h + sin * (g - h*tau);
    }

    public def getEigenValues() : Vector = eigenValuesVec;

    public def getEigenVectors() : Matrix = eigenVectorsMat; 
}
