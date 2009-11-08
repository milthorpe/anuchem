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
        if (twoE.isDirect()) { computeDirect4(twoE, density); return; }

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

    private def computeDirect(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector();
        densityOneD.make(density); // form 1D vector of density
        val tempVector:Vector{self.at(this)}  = Vector.make(noOfBasisFunctions*noOfBasisFunctions) as Vector{self.at(this)};

        val gMatrix = getMatrix();
        val temp = tempVector.getVector();

        shared var i:Int, j:Int;

        val itemsPerPlace:Int     = noOfBasisFunctions / Place.MAX_PLACES;
        val remainingItems:Int    = noOfBasisFunctions % Place.MAX_PLACES;
        val itemsOnFirstPlace:Int = noOfBasisFunctions + remainingItems;


        // TODO: x10 - parallel
        for(i=0; i<noOfBasisFunctions; i++) {
            for(j=0; j<i+1; j++) {

                tempVector.makeZero();

                val formTempVec=(startFunc:Int, endFunc:Int, i:Int, j:Int)=> {
                  var indexJ:Int, indexK1:Int, indexK2:Int;
                  var twoEIntVal1:Double, twoEIntVal2:Double, twoEIntVal3:Double;
                  var kl:Int = startFunc*noOfBasisFunctions+1;

                  for(var k:Int=startFunc; k<endFunc; k++) {
                    for(var l:Int=0; l<noOfBasisFunctions; l++) {
                        indexJ   = IntegralsUtils.ijkl2intindex(i, j, k, l);
                        indexK1  = IntegralsUtils.ijkl2intindex(i, k, j, l);
                        indexK2  = IntegralsUtils.ijkl2intindex(i, l, k, j);

                        twoEIntVal1 = twoE.compute2E(i,j,k,l);
                        if (indexJ == indexK1) twoEIntVal2 = twoEIntVal1;
                        else                   twoEIntVal2 = twoE.compute2E(i,k,j,l);

                        if (indexJ == indexK2)       twoEIntVal3 = twoEIntVal1;
                        else if (indexK1 == indexK2) twoEIntVal3 = twoEIntVal2;
                        else                         twoEIntVal3 = twoE.compute2E(i,l,k,j);

                        temp(kl++) = 2.0*twoEIntVal1 - 0.5*twoEIntVal2 - 0.5*twoEIntVal3;
                    } // end l loop
                  } // end k loop
                };
                
                shared var ii:Int = 0;
                finish for(p in Place.places) {
                    val ii_local:Int = ii;
                    val i_local:Int  = i;
                    val j_local:Int  = j;

                    if (p == Place.FIRST_PLACE)
                       async(p) formTempVec(0, itemsOnFirstPlace, i_local, j_local);
                    else
                       async(p) formTempVec(itemsOnFirstPlace + (ii_local*itemsPerPlace),
                                            itemsOnFirstPlace + (ii_local*itemsPerPlace) + itemsPerPlace,
                                            i_local, j_local);
                    ii++;
                }

                gMatrix(i, j) = gMatrix(j, i) = tempVector.dot(densityOneD);
            } // end j loop
        } // end i loop
    }

    private def computeDirect2(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        // TODO: x10 - parallel
        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;
        var twoEIntVal:Double, twoEIntVal2:Double, twoEIntValHalf:Double;
        var idx:Array[Int]{rank==1}, jdx:Array[Int]{rank==1}, 
            kdx:Array[Int]{rank==1}, ldx:Array[Int]{rank==1};
        idx = Array.make[Int]([0..8]);
        jdx = Array.make[Int]([0..8]);
        kdx = Array.make[Int]([0..8]);
        ldx = Array.make[Int]([0..8]);
        var validIdx:Array[Boolean]{rank==1} = Array.make[Boolean]([0..8]);
        validIdx(0) = true;

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
                           twoEIntVal     = twoE.compute2E(i,j,k,l);
                           twoEIntVal2    = twoEIntVal + twoEIntVal;
                           twoEIntValHalf = 0.5 * twoEIntVal;

                           setGMatrixElements(gMatrix, dMatrix, i, j, k, l,
                                               twoEIntVal2, twoEIntValHalf);

                           // special case
                           if ((i|j|k|l) == 0) continue;

                           // else this is symmetry unique integral, so need to
                           // use this value for all 8 combinations
                           // (if unique)
                           ldx(0) = l; ldx(1) = l; kdx(2) = l; kdx(3) = l;
                           idx(4) = l; idx(5) = l; jdx(6) = l; jdx(7) = l;
                           validIdx(1) = true; validIdx(2) = true;
                           validIdx(3) = true; validIdx(4) = true;
                           validIdx(5) = true; validIdx(6) = true;
                           validIdx(7) = true;

                           // filter unique elements
                           filterUniqueElements(idx, jdx, kdx, ldx, validIdx);

                           // and evaluate them
                           for(m=1; m<8; m++) {
                               if (validIdx(m)) {
                                   setGMatrixElements(gMatrix, dMatrix, 
                                               idx(m), jdx(m), kdx(m), ldx(m),
                                               twoEIntVal2, twoEIntValHalf);
                               } // end if
                           } // end for
                       } // end if
                    } // end l loop
                } // end k loop
            } // end j loop
        } // end i loop

        // half the elements
        finish ateach(var(a,b) in gMatrix.dist)
                  gMatrix(a,b) *= 0.5;
    }

    private def computeDirect3(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;
        var idx_i:Int, jdx_i:Int, kdx_i:Int, ldx_i:Int;
        var twoEIntVal:Double, twoEIntVal2:Double, twoEIntValHalf:Double;
        var idx:Array[Int]{rank==1}, jdx:Array[Int]{rank==1},
            kdx:Array[Int]{rank==1}, ldx:Array[Int]{rank==1};
        idx = Array.make[Int]([0..8]);
        jdx = Array.make[Int]([0..8]);
        kdx = Array.make[Int]([0..8]);
        ldx = Array.make[Int]([0..8]);
        var validIdx:Array[Boolean]{rank==1} = Array.make[Boolean]([0..8]);
        validIdx(0) = true;

        val molecule = twoE.getMolecule();
        val bfs = twoE.getBasisFunctions().getBasisFunctions();
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

        // TODO: x10 - parallel
        // center a
        for(a=0; a<noOfAtoms; a++) {
            aFunc = molecule.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
                iaFunc = aFunc.get(i);
                idx_i = iaFunc.getIndex();
                idx(0) = idx_i; jdx(1) = idx_i; jdx(2) = idx_i; idx(3) = idx_i;
                kdx(4) = idx_i; ldx(5) = idx_i; kdx(6) = idx_i; ldx(7) = idx_i;

                // center b
                for(b=0; b<=a; b++) {
                    bFunc = molecule.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        jbFunc = bFunc.get(j);
                        jdx_i = jbFunc.getIndex();
                        jdx(0) = jdx_i; idx(1) = jdx_i; idx(2) = jdx_i; jdx(3) = jdx_i;
                        ldx(4) = jdx_i; kdx(5) = jdx_i; ldx(6) = jdx_i; kdx(7) = jdx_i;

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            cFunc = molecule.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                kcFunc = cFunc.get(k);
                                kdx_i = kcFunc.getIndex();
                                kdx(0) = kdx_i; kdx(1) = kdx_i; ldx(2) = kdx_i; ldx(3) = kdx_i;
                                jdx(4) = kdx_i; jdx(5) = kdx_i; idx(6) = kdx_i; idx(7) = kdx_i;

                                // center d
                                for(d=0; d<=c; d++) {
                                    dFunc = molecule.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        ldFunc = dFunc.get(l);
                                        ldx_i = ldFunc.getIndex();

                                        twoEIntVal     = twoE.compute2E(iaFunc, jbFunc, kcFunc, ldFunc);
                                        twoEIntVal2    = twoEIntVal + twoEIntVal;
                                        twoEIntValHalf = 0.5 * twoEIntVal;

                                        setGMatrixElements(gMatrix, dMatrix, idx_i, jdx_i, kdx_i, ldx_i,
                                                           twoEIntVal2, twoEIntValHalf);

                                        // special case
                                        if ((idx_i|jdx_i|kdx_i|ldx_i) == 0) continue;

                                        // else this is symmetry unique integral, so need to
                                        // use this value for all 8 combinations
                                        // (if unique)
                                        ldx(0) = ldx_i; ldx(1) = ldx_i; kdx(2) = ldx_i; kdx(3) = ldx_i;
                                        idx(4) = ldx_i; idx(5) = ldx_i; jdx(6) = ldx_i; jdx(7) = ldx_i;
                                        validIdx(1) = true; validIdx(2) = true;
                                        validIdx(3) = true; validIdx(4) = true;
                                        validIdx(5) = true; validIdx(6) = true;
                                        validIdx(7) = true;

                                        // filter unique elements
                                        filterUniqueElements(idx, jdx, kdx, ldx, validIdx);

                                        // and evaluate them
                                        for(m=1; m<8; m++) {
                                          if (validIdx(m)) {
                                            setGMatrixElements(gMatrix, dMatrix,
                                                               idx(m), jdx(m), kdx(m), ldx(m),
                                                               twoEIntVal2, twoEIntValHalf);
                                          } // end if
                                        } // end m
                                    } // end l
                                } // end d
                            } // end k
                        } // end c
                    } // end j
                } // end b
            } // end i
        } // end a

        // half the elements
        finish ateach(var(x,y) in gMatrix.dist)
                  gMatrix(x,y) *= 0.5;
    }  

    private def computeDirect4(twoE:TwoElectronIntegrals{self.at(this)}, density:Density{self.at(this)}) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        // TODO: x10 - parallel
        var i:Int, j:Int, k:Int, l:Int, ij:Int, kl:Int;
        var idx:Array[Int]{rank==1}, jdx:Array[Int]{rank==1}, 
            kdx:Array[Int]{rank==1}, ldx:Array[Int]{rank==1};
        idx = Array.make[Int]([0..8]);
        jdx = Array.make[Int]([0..8]);
        kdx = Array.make[Int]([0..8]);
        ldx = Array.make[Int]([0..8]);
        
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
                         val i_loc = i, j_loc = j, k_loc = k, l_loc = l;
                         val idx_loc:Array[Int]{rank==1} = Array.make[Int]([0..8]), 
                             jdx_loc:Array[Int]{rank==1} = Array.make[Int]([0..8]), 
                             kdx_loc:Array[Int]{rank==1} = Array.make[Int]([0..8]), 
                             ldx_loc:Array[Int]{rank==1} = Array.make[Int]([0..8]);

                         for(var midx:Int=0; midx<8; midx++) {
                            idx_loc(midx) = idx(midx);
                            jdx_loc(midx) = jdx(midx);
                            kdx_loc(midx) = kdx(midx);
                            ldx_loc(midx) = ldx(midx);
                         } // end for

                         async {
                           var twoEIntVal:Double, twoEIntVal2:Double, twoEIntValHalf:Double;
                           var validIdx:Array[Boolean]{rank==1} = Array.make[Boolean]([0..8]);
                           validIdx(0) = true;

                           twoEIntVal     = twoE.compute2E(i_loc,j_loc,k_loc,l_loc);
                           twoEIntVal2    = twoEIntVal + twoEIntVal;
                           twoEIntValHalf = 0.5 * twoEIntVal;

                           setGMatrixElements(gMatrix, dMatrix, i_loc, j_loc, k_loc, l_loc,
                                               twoEIntVal2, twoEIntValHalf);

                           // special case
                           if ((i_loc|j_loc|k_loc|l_loc) != 0) {
                              // (if unique)
                              ldx_loc(0) = l_loc; ldx_loc(1) = l_loc; kdx_loc(2) = l_loc; kdx_loc(3) = l_loc;
                              idx_loc(4) = l_loc; idx_loc(5) = l_loc; jdx_loc(6) = l_loc; jdx_loc(7) = l_loc;
                              validIdx(1) = true; validIdx(2) = true;
                              validIdx(3) = true; validIdx(4) = true;
                              validIdx(5) = true; validIdx(6) = true;
                              validIdx(7) = true;
 
                              // filter unique elements
                              filterUniqueElements(idx_loc, jdx_loc, kdx_loc, ldx_loc, validIdx);

                              // and evaluate them
                              for(var m:Int=1; m<8; m++) {
                                 if (validIdx(m)) {
                                    setGMatrixElements(gMatrix, dMatrix, 
                                                idx_loc(m), jdx_loc(m), kdx_loc(m), ldx_loc(m),
                                                twoEIntVal2, twoEIntValHalf);
                                 } // end if
                              } // end for m
                           }  // end if
                         } // async block
                       } // end if
                    } // end l loop
                } // end k loop
            } // end j loop
          } // end i loop
        } // finish block

        // half the elements
        finish ateach(var(a,b) in gMatrix.dist)
                  gMatrix(a,b) *= 0.5;
    }

    /** find unique elements and mark the onces that are not */
    private def filterUniqueElements(idx:Array[Int]{rank==1}, jdx:Array[Int]{rank==1},
                                     kdx:Array[Int]{rank==1}, ldx:Array[Int]{rank==1},
                                     validIdx:Array[Boolean]{rank==1}) : void {
        var i:Int, j:Int, k:Int, l:Int, m:Int, n:Int;
        
        for(m=0; m<8; m++) {
            i = idx(m); j = jdx(m); k = kdx(m); l = ldx(m);
            for(n=m+1; n<8; n++) {
                if (i==idx(n) && j==jdx(n) && k==kdx(n) && l==ldx(n))
                    validIdx(n) = false;
            } // end for
        } // end for
    }

    /** Set the GMatrix value for a given combination */
    private def setGMatrixElements(gMatrix:Array[Double]{rank==2}, dMatrix:Array[Double]{rank==2},
                                   i:Int, j:Int, k:Int, l:Int,
                                   twoEIntVal2:Double, twoEIntValHalf:Double) : void {
        gMatrix(i,j) += dMatrix(k,l) * twoEIntVal2;
        gMatrix(k,l) += dMatrix(i,j) * twoEIntVal2;
        gMatrix(i,k) -= dMatrix(j,l) * twoEIntValHalf;
        gMatrix(i,l) -= dMatrix(j,k) * twoEIntValHalf;
        gMatrix(j,k) -= dMatrix(i,l) * twoEIntValHalf;
        gMatrix(j,l) -= dMatrix(i,k) * twoEIntValHalf;
    }
}

