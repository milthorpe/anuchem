/**
 * GMatrix.x10
 *
 * GMatrix in HF calculation
 *
 * @author: V.Ganesh
 */
 
package au.anu.edu.qm;

import x10.util.*;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

public class GMatrix extends Matrix {

    public def compute(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        if (twoE.isDirect()) { 
           val timer = new Timer(1);

           timer.start(0);
           computeDirectAsyncNew(twoE, density); 
           timer.stop(0);
           Console.OUT.println ("    Time to construct GMatrix: " + (timer.total(0) as Double) / 1e9 + " seconds");
           return; 
        } // end if 

        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector();
        densityOneD.make(density); // form 1D vector of density
        val tempVector:Vector{self.at(this)}  = Vector.make(noOfBasisFunctions*noOfBasisFunctions) as Vector{self.at(this)};

        val gMatrix = getMatrix();
        val ints = twoE.getTwoElectronIntegrals();
        val temp = tempVector.getVector();

        var i:Int, j:Int, k:Int, l:Int, 
            indexJ:Int, indexK1:Int, indexK2:Int;

        // TODO: x10 - parallel
        for(i=0; i<noOfBasisFunctions; i++) {
            for(j=0; j<i+1; j++) {

                tempVector.makeZero();
                var kl:Int = 0;

                for(k=0; k<noOfBasisFunctions; k++) {
                    for(l=0; l<noOfBasisFunctions; l++) {
                        indexJ   = IntegralsUtils.ijkl2intindex(i, j, k, l);
                        indexK1  = IntegralsUtils.ijkl2intindex(i, k, j, l);
                        indexK2  = IntegralsUtils.ijkl2intindex(i, l, k, j);
                        temp(kl++) = 2.0*ints(indexJ) - 0.5*ints(indexK1) - 0.5*ints(indexK2);
                    } // end l loop
                } // end k loop

                gMatrix(i, j) = gMatrix(j, i) = tempVector.dot(densityOneD);
            } // end j loop
        } // end i loop  
    }

    private def computeDirectSerialOld(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = Matrix.make(N) as Matrix{self.at(this)};
        val kMat = Matrix.make(N) as Matrix{self.at(this)};

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;
        var idx_i:Int, jdx_i:Int, kdx_i:Int, ldx_i:Int;
        val idx = Rail.make[Int](8);
        val jdx = Rail.make[Int](8);
        val kdx = Rail.make[Int](8);
        val ldx = Rail.make[Int](8);

        for(i=0; i<N; i++) {
            idx(0) = i; jdx(1) = i; jdx(2) = i; idx(3) = i;
            kdx(4) = i; ldx(5) = i; kdx(6) = i; ldx(7) = i;
            for(j=0; j<(i+1); j++) {
                ij = i * (i+1) / 2+j;
                jdx(0) = j; idx(1) = j; idx(2) = j; jdx(3) = j;
                ldx(4) = j; kdx(5) = j; ldx(6) = j; kdx(7) = j;
                for(k=0; k<N; k++) {
                    kdx(0) = k; kdx(1) = k; ldx(2) = k; ldx(3) = k;
                    jdx(4) = k; jdx(5) = k; idx(6) = k; idx(7) = k;
                    for(l=0; l<(k+1); l++) {
                        kl = k * (k+1) / 2+l;
                        if (ij >= kl) { 
                           val twoEIntVal = twoE.compute2E(i,j,k,l);

                           setJKMatrixElements(jMatrix, kMatrix, dMatrix, i, j, k, l, twoEIntVal);

                           // special case
                           if ((i|j|k|l) == 0) continue;

                           // else this is symmetry unique integral, so need to
                           // use this value for all 8 combinations
                           // (if unique)
                           ldx(0) = l; ldx(1) = l; kdx(2) = l; kdx(3) = l;
                           idx(4) = l; idx(5) = l; jdx(6) = l; jdx(7) = l;
                           val validIdx = Rail.make[Boolean](8, (Int)=>true);

                           // filter unique elements
                           filterUniqueElements(idx, jdx, kdx, ldx, validIdx);

                           // and evaluate them
                           for(m=1; m<8; m++) {
                               if (validIdx(m)) {
                                   setJKMatrixElements(jMatrix, kMatrix, dMatrix, 
                                                       idx(m), jdx(m), kdx(m), ldx(m),
                                                       twoEIntVal);
                               } // end if
                           } // end for
                       } // end if                        
                    } // end l loop
                } // end k loop
            } // end j loop
        } // end i loop

        // form the G matrix
        finish ateach(val(x,y) in gMatrix.dist)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));

    }

    private def computeDirectAsyncOld(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = Matrix.make(N) as Matrix{self.at(this)};
        val kMat = Matrix.make(N) as Matrix{self.at(this)};

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;
        val idx = Rail.make[Int](8);
        val jdx = Rail.make[Int](8);
        val kdx = Rail.make[Int](8);
        val ldx = Rail.make[Int](8);

        finish {
          for(i=0; i<N; i++) {
            idx(0) = i; jdx(1) = i; jdx(2) = i; idx(3) = i;
            kdx(4) = i; ldx(5) = i; kdx(6) = i; ldx(7) = i;
            for(j=0; j<(i+1); j++) {
                ij = i * (i+1) / 2+j;
                jdx(0) = j; idx(1) = j; idx(2) = j; jdx(3) = j;
                ldx(4) = j; kdx(5) = j; ldx(6) = j; kdx(7) = j;
                for(k=0; k<N; k++) {
                    kdx(0) = k; kdx(1) = k; ldx(2) = k; ldx(3) = k;
                    jdx(4) = k; jdx(5) = k; idx(6) = k; idx(7) = k;
                    for(l=0; l<(k+1); l++) {
                        kl = k * (k+1) / 2+l;
                        if (ij >= kl) { 
			   val i_loc = i, j_loc = j , k_loc = k, l_loc = l;
                           val idx_loc = Rail.make[Int](8);
                           val jdx_loc = Rail.make[Int](8);
                           val kdx_loc = Rail.make[Int](8);
                           val ldx_loc = Rail.make[Int](8);

                           for(var midx:Int=0; midx<8; midx++) {
                               idx_loc(midx) = idx(midx);
                               jdx_loc(midx) = jdx(midx);
                               kdx_loc(midx) = kdx(midx);
                               ldx_loc(midx) = ldx(midx);
                           } // end for

                           async {
                                val twoEIntVal = twoE.compute2E(i_loc, j_loc, k_loc, l_loc);
                                val validIdx = Rail.make[Boolean](8, (Int)=>true);

                                setJKMatrixElements(jMatrix, kMatrix, dMatrix, i_loc, j_loc, k_loc, l_loc, twoEIntVal);

                                // if not special case
                                if ((i_loc|j_loc|k_loc|l_loc) != 0) {
                                    // else this is symmetry unique integral, so need to
                                    // use this value for all 8 combinations
                                    // (if unique)
                                    ldx_loc(0) = l_loc; ldx_loc(1) = l_loc; kdx_loc(2) = l_loc; kdx_loc(3) = l_loc;
                                    idx_loc(4) = l_loc; idx_loc(5) = l_loc; jdx_loc(6) = l_loc; jdx_loc(7) = l_loc;

                                    // filter unique elements
                                    filterUniqueElements(idx_loc, jdx_loc, kdx_loc, ldx_loc, validIdx);

                                    // and evaluate them
                                    for(var m:Int=1; m<8; m++) {
                                        if (validIdx(m)) {
                                            setJKMatrixElements(jMatrix, kMatrix, dMatrix,
                                                                idx_loc(m), jdx_loc(m), kdx_loc(m), ldx_loc(m),
                                                                twoEIntVal);
                                        } // end if
                                    } // end m
                                } // end if
                           } // async block
                       } // end if                        
                    } // end l loop
                } // end k loop
            } // end j loop
          } // end i loop
        } // finish

        // form the G matrix
        finish ateach(val(x,y) in gMatrix.dist)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));
    }

    private def computeDirectSerialNew(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = Matrix.make(N) as Matrix{self.at(this)};
        val kMat = Matrix.make(N) as Matrix{self.at(this)};

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        var i:Int, j:Int, k:Int, l:Int;

        val molecule = twoE.getMolecule();
        val bfs = twoE.getBasisFunctions().getBasisFunctions();
        val shellList = twoE.getBasisFunctions().getShellList();
        val noOfBasisFunctions = bfs.size();

        val noOfAtoms = molecule.getNumberOfAtoms();
        var a:Int, b:Int, c:Int, d:Int;
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
            kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};

        // center a
        for(a=0; a<noOfAtoms; a++) {
            aFunc = molecule.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
                iaFunc = aFunc.get(i);

                // center b
                for(b=0; b<=a; b++) {
                    bFunc = molecule.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        jbFunc = bFunc.get(j);

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            cFunc = molecule.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                kcFunc = cFunc.get(k);

                                // center d
                                for(d=0; d<=c; d++) {
                                    dFunc = molecule.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        ldFunc = dFunc.get(l);

            	                        twoE.compute2EAndRecord(iaFunc, jbFunc, kcFunc, ldFunc, 
                                                                shellList, jMat, kMat, density);
				    } // end l
				} // center d
                            } // end k
                        } // center c
                    } // end j
                } // center b
            } // end i
        } // center a
     
        // form the G matrix
        finish ateach(val(x,y) in gMatrix.dist)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    private def computeDirectAsyncNew(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = Matrix.make(N) as Matrix{self.at(this)};
        val kMat = Matrix.make(N) as Matrix{self.at(this)};

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        var i:Int, j:Int, k:Int, l:Int;

        val molecule = twoE.getMolecule();
        val bfs = twoE.getBasisFunctions().getBasisFunctions();
        val shellList = twoE.getBasisFunctions().getShellList();
        val noOfBasisFunctions = bfs.size();

        val noOfAtoms = molecule.getNumberOfAtoms();
        var a:Int, b:Int, c:Int, d:Int;
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
            kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};

        finish {
          // center a
          for(a=0; a<noOfAtoms; a++) {
            aFunc = molecule.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
                iaFunc = aFunc.get(i);

                // center b
                for(b=0; b<=a; b++) {
                    bFunc = molecule.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        jbFunc = bFunc.get(j);

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            cFunc = molecule.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                kcFunc = cFunc.get(k);

                                // center d
                                for(d=0; d<=c; d++) {
                                    dFunc = molecule.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        ldFunc = dFunc.get(l);

                                        val iaFunc_loc = iaFunc;
                                        val jbFunc_loc = jbFunc;
                                        val kcFunc_loc = kcFunc;
                                        val ldFunc_loc = ldFunc;

                                        async {                                         
              	                           twoE.compute2EAndRecord(iaFunc_loc, jbFunc_loc, kcFunc_loc, ldFunc_loc, 
                                                                   shellList, jMat, kMat, density);
                                        }
				    } // end l
				} // center d
                            } // end k
                        } // center c
                    } // end j
                } // center b
            } // end i
          } // center a
        } // finish
      
        // form the G matrix
        finish ateach(val(x,y) in gMatrix.dist)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    /** find unique elements and mark the onces that are not */
    private def filterUniqueElements(idx:Rail[Int]!, jdx:Rail[Int]!,
                                     kdx:Rail[Int]!, ldx:Rail[Int]!,
                                     validIdx:Rail[Boolean]!) : void {
        var i:Int, j:Int, k:Int, l:Int, m:Int, n:Int;
        
        for(m=0; m<8; m++) {
            i = idx(m); j = jdx(m); k = kdx(m); l = ldx(m);
            for(n=m+1; n<8; n++) {
                if (i==idx(n) && j==jdx(n) && k==kdx(n) && l==ldx(n))
                    validIdx(n) = false;
            } // end for
        } // end for
    }

    /** Set the J and K value for a given combination */
    private def setJKMatrixElements(jMatrix:Array[Double]{rank==2}, kMatrix:Array[Double]{rank==2},
                                    dMatrix:Array[Double]{rank==2},
                                    i:Int, j:Int, k:Int, l:Int, twoEIntVal:Double) : void {
        val v1 = dMatrix(k,l) * twoEIntVal;
        val v2 = dMatrix(i,j) * twoEIntVal;
        val v3 = dMatrix(j,l) * twoEIntVal;
        val v4 = dMatrix(j,k) * twoEIntVal;
        val v5 = dMatrix(i,l) * twoEIntVal;
        val v6 = dMatrix(i,k) * twoEIntVal;

        atomic {
          jMatrix(i,j) += v1;
          jMatrix(k,l) += v2;
          kMatrix(i,k) += v3;
          kMatrix(i,l) += v4;
          kMatrix(j,k) += v5;
          kMatrix(j,l) += v6;
        } // atomic
    }
}

