package au.anu.edu.qm;

import x10x.matrix.Matrix;

public class Fock extends Matrix {
    public def this(siz:Int) {
       super(siz);
    }

    public def compute(hCore:HCore, gMatrix:GMatrix) : void {
        mat = hCore.add(gMatrix).getMatrix();
    }
}

