package au.edu.anu.qm;

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;

import x10.util.ArrayList;
import x10.util.ArrayUtils;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;

public class TestDGEMM{
    public static def main(args:Array[String](1)){
        Console.OUT.println("Hello World!");
        val timer = au.edu.anu.util.new StatisticalTimer(4);

        val N=700; val nOrbital=25; val roK=196; // if numbers are big => GC Warning: Out of Memory!  Returning NIL!, if it's really big ==> Segmentation fault (core dumped)
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

        DenseMatrixBLAS.compMultTrans(mos, auxIntMat, halfAuxMat, [nOrbital,N*roK, N], false);
        val halfAuxMat2 = new DenseMatrix(roK*nOrbital, N, halfAuxMat.d);
        DenseMatrixBLAS.compTransMult(halfAuxMat2, halfAuxMat2, kMatrix, [N, N, roK*nOrbital], true);
        Console.OUT.println("DGEMM done!");
    }
}
