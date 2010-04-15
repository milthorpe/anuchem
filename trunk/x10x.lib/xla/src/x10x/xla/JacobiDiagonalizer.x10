package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

/**
 * Based on JacobiDiagonalizer.java (org.meta.math.la.JacobiDiagonalizer)
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public class JacobiDiagonalizer {
    var eigenValuesVec:Vector!; 
    var eigenVectorsMat:Matrix!;
    var eigenValues:Array[Double]{rank==1};
    var eigenVectors:Array[Double]{rank==2};
    var cos:Double, sin:Double, tau:Double;

    global val maxIterations : Int = 100;

    public def this() { }

    public def diagonalize(mat:Matrix!) : void {
       val matrix = mat.getMatrix();
       val n  = mat.getRowCount();
       val n1 = n-1;

       eigenVectorsMat = new Matrix(Dist.make([0..n1,0..n1])) as Matrix!;
       eigenVectorsMat.makeIdentity();

       eigenVectors = eigenVectorsMat.getMatrix();

       // clone the matrix (make it single place), do not tamper the actual matrix
       val aMat = new Matrix(Dist.make([0..n1,0..n1])) as Matrix!;
       val a = aMat.getMatrix();

       finish foreach(val(i,j) in matrix.region) { 
                 a(i, j) = matrix(i, j); 
       }

       eigenValuesVec = new Vector(Dist.make([0..n1])) as Vector!; 
       eigenValues = eigenValuesVec.getVector();
      
       val bVec = new Vector(Dist.make([0..n1])) as Vector!;
       val zVec = new Vector(Dist.make([0..n1])) as Vector!;
       val b = bVec.getVector();
       val z = zVec.getVector();

       // do diagonalization at Place.FIRST_PLACE
       finish foreach(val(i,j) in a.region) {
           if (i == j) { 
               eigenValues(i) = b(i) = a(i,i);
               z(i) = 0.0;
           } // end if
       } // end ateach

       var zeroTolerance:Double = 0.0;

       for(var sweepsIdx:Int = 0; sweepsIdx<maxIterations; sweepsIdx++) {
          // sum off diagonal elements of A
          val sum = aMat.sumOffDiagonal();
          
          if (sum == 0.0) break;  // if off diagonal elements are zero, we stop

          if (sweepsIdx < 3) zeroTolerance = 0.2 * sum / (n*n);
          else               zeroTolerance = 0.0;

          val sweeps = sweepsIdx;
          val zeroTol = zeroTolerance;

          for(var ip:Int = 0; ip<n-1; ip++) { for(var iq:Int = ip+1; iq<n; iq++) {

          // finish foreach(plc in a.dist.places()) { at(plc) { for(val(ip, iq) in a.dist.get(plc)) {
             // x10.io.Console.OUT.println(ip + " , " + iq);
             // x10.io.Console.OUT.println(a.dist.contains(Point.make(ip, iq)));

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

                     cos = 1.0 / Math.sqrt(1.0 + t*t);
                     sin = t * cos;
                     tau = sin / (1.0 + cos);
                     h   = t * a(ip, iq);

                     z(ip) -= h; z(iq) += h;
                     eigenValues(ip) -= h; eigenValues(iq) += h;

                     a(ip,iq) = 0.0;

                     /**
                     for(val(i,j) in a.dist) if (i==j && j<ip)             
                         doRotate(a, j, ip, j, iq, sin, tau);
                     for(val(i,j) in a.dist) if (i==j && j>=ip+1 && j<iq)  
                         doRotate(a, ip, j, j, iq, sin, tau);
                     for(val(i,j) in a.dist) if (i==j && j>=iq+1 && j<n)   
                         doRotate(a, ip, j, iq, j, sin, tau);
                     for(val(i,j) in eigenVectors.dist) if (i==j)          
                         doRotate(eigenVectors, j, ip, j, iq, sin, tau);
                     **/

                     for(var j:Int = 0; j<ip; j++)    doRotate(a, j, ip, j, iq, sin, tau);
                     for(var j:Int = ip+1; j<iq; j++) doRotate(a, ip, j, j, iq, sin, tau);
                     for(var j:Int = iq+1; j<n; j++)  doRotate(a, ip, j, iq, j, sin, tau);
                     for(var j:Int = 0; j<n; j++)     doRotate(eigenVectors, j, ip, j, iq, sin, tau);
                   } // end if
                } // end if
             } // end if
           // }}} // end foreach
           }} // end for

          finish foreach(val(ip) in b.region) {
             b(ip) += z(ip);
             eigenValues(ip) = b(ip);
             z(ip) = 0.0;
          } // end ateach
       } // end for

       sortEigenValues(); 
       eigenVectorsMat = eigenVectorsMat.transpose();
    }

    private static def doRotate(a:Array[Double]{rank==2}, i:Int, j:Int, k:Int, l:Int,
                                sin:Double, tau:Double) : void {
       var g:Double = a(i,j);
       var h:Double = a(k,l);

       a(i,j) = g - sin * (h + g*tau);
       a(k,l) = h + sin * (g - h*tau);
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

    public def getEigenValues() : Vector! = eigenValuesVec;

    public def getEigenVectors() : Matrix! = eigenVectorsMat; 
}

