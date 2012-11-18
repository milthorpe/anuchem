package au.edu.anu.qm;

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;

import x10.util.ArrayList;
import x10.util.ArrayUtils;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;

import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;

public class TestDGEMM{
    public static def main(args:Array[String](1)){
        Console.OUT.println("Hello World!");
        val timer = new StatisticalTimer(4);

        val N=1400; val nOrbital=25; val roK=196; // if numbers are big => GC Warning: Out of Memory!  Returning NIL!, if it's really big ==> Segmentation fault (core dumped)
        val kMatrix=new DenseMatrix(N,N);
        val auxIntMat=new DenseMatrix(N*roK,N);
        val halfAuxMat=new DenseMatrix(nOrbital,N*roK);
        val mos=new DenseMatrix(nOrbital,N);
        Console.OUT.println("Variables allocated.");      

        for ([i,j] in (0..(N*roK-1))*(0..(N-1)))
            auxIntMat(i,j)=i*123+j;        
        for ([i,j] in (0..(nOrbital-1))*(0..(N*roK-1)))
            halfAuxMat(i,j)=i*123+j;
        for ([i,j] in (0..(nOrbital-1))*(0..(N-1)))
            mos(i,j)=i*123+j;
        Console.OUT.println("Variables assigned.");

        timer.start(0);  
        DenseMatrixBLAS.compMultTrans(mos, auxIntMat, halfAuxMat, [nOrbital,N*roK, N], false);
        timer.stop(0);
        Console.OUT.println("DGEMM done! " + (timer.last(0) as Double)/1e9 + "s");
 
        timer.start(1);  
        val halfAuxMat2 = new DenseMatrix(roK*nOrbital, N, halfAuxMat.d);
        timer.stop(1);
        Console.OUT.println("Cast done! " + (timer.last(1) as Double)/1e9 + "s");

        timer.start(2);  
        DenseMatrixBLAS.compTransMult(halfAuxMat2, halfAuxMat2, kMatrix, [N, N, roK*nOrbital], true);
        timer.stop(2); 
        Console.OUT.println("DGEMM done! "+ (timer.last(2) as Double)/1e9 + "s");
    }
}
