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
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;

public class GMatrix extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    /** top level method to form the G Matrix, depending on gMatType appropriate functions are called */
    public def compute(twoE:TwoElectronIntegrals!, density:Density!, gMatType:Int) : void {
        if (twoE.isDirect()) { 
           val timer = new Timer(1);

           timer.start(0);
           switch(gMatType) {
           case 0:
               Console.OUT.println("   GMatrix.computeDirectSerialOld: " + gMatType);
               computeDirectSerialOld(twoE, density); 
               break;
           case 1:
               Console.OUT.println("   GMatrix.computeDirectAsyncOld: " + gMatType);
               computeDirectAsyncOld(twoE, density); 
               break;
           case 2:
               Console.OUT.println("   GMatrix.computeDirectAsyncOldNoAtomic: " + gMatType);
               computeDirectAsyncOldNoAtomic(twoE, density); 
               break;
           case 3:
               Console.OUT.println("   GMatrix.computeDirectSerialNew: " + gMatType);
               computeDirectSerialNew(twoE, density); 
               break;
           case 4:
               Console.OUT.println("   GMatrix.computeDirectLowMemNewNoAtomic: " + gMatType);
               computeDirectLowMemNewNoAtomic(twoE, density); 
               break;
	   case 5:
               Console.OUT.println("   GMatrix.computeDirectOldMultiPlaceNoAtomic: " + gMatType);
 	       computeDirectOldMultiPlaceNoAtomic(twoE, density);
	       break;
           case 6:
               Console.OUT.println("   GMatrix.computeDirectNewMultiPlaceNoAtomic: " + gMatType);
               computeDirectNewMultiPlaceNoAtomic(twoE, density);
               break;
           case 7:
               Console.OUT.println("   GMatrix.computeDirectOldMultiPlaceTaskPool: " + gMatType);
               computeDirectOldMultiPlaceTaskPool(twoE, density);
               break;
           case 8:
               Console.OUT.println("   GMatrix.computeDirectOldMultiPlaceStatic: " + gMatType);
               computeDirectOldMultiPlaceStatic(twoE, density);
               break;
           case 9:
               Console.OUT.println("   GMatrix.computeDirectNewMultiPlaceStatic: " + gMatType);
               computeDirectNewMultiPlaceStatic(twoE, density);
               break;
           case 10:
               Console.OUT.println("   GMatrix.computeDirectMultiPlaceNewFuture: " + gMatType);
               computeDirectMultiPlaceNewFuture(twoE, density);
               break;
           default:
               Console.OUT.println("   GMatrix.computeDirectSerialOld: " + gMatType);
               computeDirectSerialOld(twoE, density); 
               break;
           } // end switch .. case
           timer.stop(0);
           Console.OUT.println ("    Time to construct GMatrix: " + (timer.total(0) as Double) / 1e9 + " seconds");
           return; 
        } // end if 

        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector(density); // form 1D vector of density
        val tempVector  = new Vector(noOfBasisFunctions*noOfBasisFunctions) as Vector!;

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

    /* Plain serial version */
    private def computeDirectSerialOld(twoE:TwoElectronIntegrals!, density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = new Matrix(N) as Matrix!;
        val kMat = new Matrix(N) as Matrix!;

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
        finish foreach(val(x,y) in gMatrix.region)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));

    }

    /* Simple, but crashes for larger test case, and uses atomic */
    private def computeDirectAsyncOld(twoE:TwoElectronIntegrals!, density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = new Matrix(N) as Matrix!;
        val kMat = new Matrix(N) as Matrix!;

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
        finish foreach(val(x,y) in gMatrix.region)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));
    }

    /* No atomic used, local GMatrix per thread, but introduces artificial sync */
    private def computeDirectAsyncOldNoAtomic(twoE:TwoElectronIntegrals!, 
                                              density:Density!) : void {
        val N = density.getRowCount();
         
        val molecule = twoE.getMolecule();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;

        val idx = Rail.make[Int](8);
        val jdx = Rail.make[Int](8);
        val kdx = Rail.make[Int](8);
        val ldx = Rail.make[Int](8);

        val computeInst = Rail.make[ComputeOld!](Runtime.INIT_THREADS);

        for(i=0; i<Runtime.INIT_THREADS; i++) {
            computeInst(i) = new ComputeOld(new TwoElectronIntegrals(twoE.getBasisFunctions(), molecule, true), density);
        } // end for

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
                           
                           var setIt:Boolean = false;

                           outer: while(!setIt) {
                               for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                                   setIt = computeInst(ix).setValue(i, j, k, l, idx, jdx, kdx, ldx);
                                   if (setIt) {
                                      val ix_loc = ix;
                                      async computeInst(ix_loc).compute();
                                      break outer;
                                   } // end if
                               } // end for
                           } // end while
                       } // end if                        
                    } // end l loop
                } // end k loop
            } // end j loop
          } // end i loop
        } // finish

        // form the G matrix
        for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
             val jMatrix = computeInst(ix).getJMat().getMatrix();
             val kMatrix = computeInst(ix).getKMat().getMatrix();

             finish foreach(val(x,y) in gMatrix.region) {
                   gMatrix(x,y) += jMatrix(x,y) - (0.25*kMatrix(x,y));
             } // finish
        } // end for
    }

    /* multi place version of the above code */
    private def computeDirectOldMultiPlaceNoAtomic(twoE:TwoElectronIntegrals!, 
                                                   density:Density!) : void {
        val N = density.getRowCount();
         
        val molecule = twoE.getMolecule();
        val basisFunctions = twoE.getBasisFunctions();

        makeZero();

        val gMatrix = getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceOld](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);

        val timer = new Timer(3);

        timer.start(0);
        i = 0;
        for(place in Place.places) {
           computeInst(i) = at(place) { 
                return new ComputePlaceOld(molecule, basisFunctions, density, place); 
           };
           i++;
        }
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        finish {
          for(i=0; i<N; i++) {
            for(j=0; j<(i+1); j++) {
                ij = i * (i+1) / 2+j;
                for(k=0; k<N; k++) {
                    for(l=0; l<(k+1); l++) {
                        kl = k * (k+1) / 2+l;
                        if (ij >= kl) { 
                           
                           var setIt:Boolean = false;
                           
                           val i_l = i, j_l = j, k_l = k, l_l = l;

                           outer: while(!setIt) {
                               for(comp_loc in computeInst) {
                                   setIt = at(comp_loc) { comp_loc.setValue(i_l, j_l, k_l, l_l) };

                                   if (setIt) break outer;
                               } // end for
                           } // end while
                       } // end if                        
                    } // end l loop
                } // end k loop
            } // end j loop
          } // end i loop
        } // finish
        timer.stop(1);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(1) as Double) / 1e9 + " seconds"); 


        timer.start(2);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };
             
                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for
        
        timer.stop(2);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(2) as Double) / 1e9 + " seconds"); 
    }

    /** multi place, but uses a static load balance */
    private def computeDirectOldMultiPlaceStatic(twoE:TwoElectronIntegrals!, 
                                                 density:Density!) : void {
        val N = density.getRowCount();
         
        val molecule = twoE.getMolecule();
        val basisFunctions = twoE.getBasisFunctions();

        makeZero();

        val gMatrix = getMatrix();

        var i:Int, j:Int;

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceOldDirect](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);

        val timer = new Timer(3);

        timer.start(0);
        i = 0;
        for(place in Place.places) {
           computeInst(i) = at(place) { return new ComputePlaceOldDirect(molecule, basisFunctions, density, place); };
           i++;
        }

        val workPerPlace = Rail.make[Int](nPlaces, (Int)=>0);
        
        i = N;
        while(i > 0) {
           for(j=0; j<nPlaces; j++) {
              workPerPlace(j) += (i <= 0) ? 0 : 1;
              i--;
           } // end for
        } // end while
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        Console.OUT.print("\tWorks units per place: ");
        for(j=0; j<nPlaces; j++) Console.OUT.print(workPerPlace(j) + ", ");
        Console.OUT.println(" "); 

        timer.start(1);
        finish {
          
          var curInd:Int = 0;
          i = 0;
          for (place in Place.places) {
              val comp_loc = computeInst(i);
              val st = curInd;
              val ed = st + workPerPlace(i);

              async at(comp_loc) { comp_loc.compute(st, ed); };
              curInd = ed;
              i++;
          } // end for
        } // finish
        timer.stop(1);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(1) as Double) / 1e9 + " seconds"); 

        timer.start(2);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };
             
                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for
        
        timer.stop(2);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(2) as Double) / 1e9 + " seconds"); 
    }

    /* multi place with tak pool interface */
    private def computeDirectOldMultiPlaceTaskPool(twoE:TwoElectronIntegrals!, 
                                                   density:Density!) : void {
        val N = density.getRowCount();
         
        val molecule = twoE.getMolecule();
        val basisFunctions = twoE.getBasisFunctions();

        makeZero();

        val gMatrix = getMatrix();

        var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceOldPool](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);

        val timer = new Timer(4);

        timer.start(0);
        i = 0;
        for(place in Place.places) {
           computeInst(i) = at(place) { 
                return new ComputePlaceOldPool(molecule, basisFunctions, density, place); 
           };
           i++;
        }
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        finish {
          ////
          // for(comp_loc in computeInst) async at(comp_loc) { comp_loc.start(); };
          ////

          timer.start(2);
          var ix:Int = 0;
          for(i=0; i<N; i++) {
            for(j=0; j<(i+1); j++) {
                ij = i * (i+1) / 2+j;
                for(k=0; k<N; k++) {
                    for(l=0; l<(k+1); l++) {
                        kl = k * (k+1) / 2+l;
                        if (ij >= kl) { 
                           val i_l = i, j_l = j, k_l = k, l_l = l;

                           if (ix >= nPlaces) ix = 0;

                           val comp_loc = computeInst(ix);

                           at(comp_loc) { comp_loc.pushBlockIdx(BlockIndices(i_l, j_l, k_l, l_l)); };

                           ix++;
                       } // end if                        
                    } // end l loop
                } // end k loop
            } // end j loop
          } // end i loop
          timer.stop(2);
          Console.OUT.println("\tTime for setting up the queue: " + (timer.total(2) as Double) / 1e9 + " seconds"); 
          
          ////
          for(comp_loc in computeInst) async at(comp_loc) { comp_loc.start(); };
          for(comp_loc in computeInst) async at(comp_loc) { comp_loc.noMoreTasks = true; };
          ////
        } // finish
        timer.stop(1);
        Console.OUT.println("\tTotal time for computation: " + (timer.total(1) as Double) / 1e9 + " seconds"); 


        timer.start(3);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };
             
                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for
        
        timer.stop(3);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(2) as Double) / 1e9 + " seconds"); 
    }

    private def computeDirectSerialNew(twoE:TwoElectronIntegrals!, density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = new Matrix(N) as Matrix!;
        val kMat = new Matrix(N) as Matrix!;

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
        finish foreach(val(x,y) in gMatrix.region)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    private def computeDirectLowMemNewNoAtomic(twoE:TwoElectronIntegrals!, 
                                               density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val molecule = twoE.getMolecule();
        val bfs = twoE.getBasisFunctions().getBasisFunctions();
        val shellList = twoE.getBasisFunctions().getShellList();
        val noOfBasisFunctions = bfs.size();

        val noOfAtoms = molecule.getNumberOfAtoms();

        // create another future to feed in the futures created above
        var i:Int, j:Int, k:Int, l:Int;
        var a:Int, b:Int, c:Int, d:Int;
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
              kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

        val computeInst = Rail.make[ComputeNew!](Runtime.INIT_THREADS);

        for(i=0; i<Runtime.INIT_THREADS; i++) {
            computeInst(i) = new ComputeNew(new TwoElectronIntegrals(twoE.getBasisFunctions(), molecule, true),
                                            shellList, density);
        } // end for

        // center a
        finish {
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

                                       var setIt:Boolean = false;
                                       
                                       outer: while(!setIt) {
                                           for(var idx:Int=0; idx<Runtime.INIT_THREADS; idx++) { 
                                               setIt = computeInst(idx).setValue(iaFunc, jbFunc, kcFunc, ldFunc);                                               
                                               if (setIt) {
                                                  val idx_loc = idx;
                                                  async computeInst(idx_loc).compute();
                                                  break outer;
                                               } // end if
                                           } // end for
                                       } // end while
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
        for(var idx:Int=0; idx<Runtime.INIT_THREADS; idx++) { 
             val jMatrix = computeInst(idx).getJMat().getMatrix();
             val kMatrix = computeInst(idx).getKMat().getMatrix();
             
             finish foreach(val(x,y) in gMatrix.region) {
                   gMatrix(x,y) += jMatrix(x,y) - (0.25*kMatrix(x,y));     
             } // finish
        } // end for
    }

    private def computeDirectNewMultiPlaceNoAtomic(twoE:TwoElectronIntegrals!, 
                                                   density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val molecule = twoE.getMolecule();
        val basisFunctions = twoE.getBasisFunctions();
        val bfs = twoE.getBasisFunctions().getBasisFunctions();
        val shellList = twoE.getBasisFunctions().getShellList();
        val noOfBasisFunctions = bfs.size();

        val noOfAtoms = molecule.getNumberOfAtoms();

        // create another future to feed in the futures created above
        var i:Int, j:Int, k:Int, l:Int;
        var a:Int, b:Int, c:Int, d:Int;
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
              kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
              dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceNew](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);

        val timer = new Timer(3);

        timer.start(0);
        i = 0;
        for(place in Place.places) {
           computeInst(i) = at(place) { 
                return new ComputePlaceNew(molecule, basisFunctions, density, place); 
           };
           i++;
        }
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        // center a
        finish {
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

				       var setIt:Boolean = false;
                           
                                       val a_l = a, b_l = b, c_l = c, d_l = d;
                                       val i_l = i, j_l = j, k_l = k, l_l = l;

                                       outer: while(!setIt) {
                                         for(comp_loc in computeInst) {
                                           setIt = at(comp_loc) { comp_loc.setValue(a_l, b_l, c_l, d_l, i_l, j_l, k_l, l_l) };

                                           if (setIt) break outer;
                                         } // end for
                                       } // end while                               
				   } // end l
			       } // center d
                           } // end k
                       } // center c
                   } // end j
               } // center b
            } // end i
          } // center a
        } // finish
        timer.stop(1);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(1) as Double) / 1e9 + " seconds"); 

        timer.start(2);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };
             
                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for
        
        timer.stop(2);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(0) as Double) / 1e9 + " seconds"); 
    }

    private def computeDirectNewMultiPlaceStatic(twoE:TwoElectronIntegrals!, 
                                                 density:Density!) : void {
        val N = density.getRowCount();
         
        val molecule = twoE.getMolecule();
        val basisFunctions = twoE.getBasisFunctions();

        makeZero();

        val gMatrix = getMatrix();

        var i:Int, j:Int;

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceNewDirect](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);

        val timer = new Timer(3);

        timer.start(0);
        i = 0;
        for(place in Place.places) {
           computeInst(i) = at(place) { return new ComputePlaceNewDirect(molecule, basisFunctions, density, place); };
           i++;
        }

        val workPerPlace = Rail.make[Int](nPlaces, (Int)=>0);
        
        i = molecule.getNumberOfAtoms();
        while(i > 0) {
           for(j=0; j<nPlaces; j++) {
              workPerPlace(j) += (i <= 0) ? 0 : 1;
              i--;
           } // end for
        } // end while
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        Console.OUT.print("\tWorks units per place: ");
        for(j=0; j<nPlaces; j++) Console.OUT.print(workPerPlace(j) + ", ");
        Console.OUT.println(" "); 

        timer.start(1);
        finish {
          
          var curInd:Int = 0;
          i = 0;
          for (place in Place.places) {
              val comp_loc = computeInst(i);
              val st = curInd;
              val ed = st + workPerPlace(i);

              async at(comp_loc) { comp_loc.compute(st, ed); };
              curInd = ed;
              i++;
          } // end for
        } // finish
        timer.stop(1);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(1) as Double) / 1e9 + " seconds"); 

        timer.start(2);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };
             
                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for
        
        timer.stop(2);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(2) as Double) / 1e9 + " seconds"); 
    }


    // TODO: need to use shared variable instead
    private global val G = Rail.make[Int](1, (Int)=>0);
    private global val PIdx = Rail.make[Int](1, (Int)=>0);

    /** Code snippet 3, Bernholdt paper  */
    private def computeDirectMultiPlaceNewFuture(twoE:TwoElectronIntegrals!, density:Density!) : void {
        // init counters
        G(0) = 0; PIdx(0) = 0;

        // init other variables
        val N = density.getRowCount();

        makeZero();

        val gMatrix = getMatrix();
        val dMatrix = density.getMatrix();

        val jMat = new Matrix(N) as Matrix!;
        val kMat = new Matrix(N) as Matrix!;

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        val molecule = twoE.getMolecule();
        val bfs = twoE.getBasisFunctions();
        val shellList = twoE.getBasisFunctions().getShellList();
        val noOfBasisFunctions = bfs.getBasisFunctions().size();

        val noOfAtoms = molecule.getNumberOfAtoms();
        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceNewFuture](nPlaces);

        Console.OUT.println("\tNo. of places: " + nPlaces);
        Console.OUT.println("\tNo. of threads per place: " + Runtime.INIT_THREADS);
  
        val timer = new Timer(2);

        timer.start(0);
    
        // center a
        finish foreach(place in Place.places) {
          async at(place) {
            var myG:Int = 0;
            var L:Int = 0;

            var i:Int, j:Int, k:Int, l:Int;
            var a:Int, b:Int, c:Int, d:Int;
            var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;

            Console.OUT.println("Copying data..." + here);
            // TODO: better way to pass global data instead?
            // make local copies of data 
            val mol_loc = new Molecule[QMAtom]() as Molecule[QMAtom]!;

            for(i=0; i<noOfAtoms; i++) {
                val i_loc = i;
                val sym  = at(molecule) { molecule.getAtom(i_loc).symbol };
                val x    = at(molecule) { molecule.getAtom(i_loc).centre.i };
                val y    = at(molecule) { molecule.getAtom(i_loc).centre.j };
                val z    = at(molecule) { molecule.getAtom(i_loc).centre.k };

                mol_loc.addAtom(new QMAtom(sym, new Point3d(x, y, z)));
            } // end for

            val basisName = at(bfs) { bfs.getBasisName() };
            val bas_loc = new BasisFunctions(mol_loc, basisName, "basis");

            val den_loc = new Density(N);
            val den_loc_mat = den_loc.getMatrix();
            val den_rem = at(density) { density.getValRail() };
            var ii:Int;

            ii = 0 ;
            for(i=0; i<N; i++)
              for(j=0; j<N; j++)
                  den_loc_mat(i, j) = den_rem(ii++);

            val twoE_loc = new TwoElectronIntegrals(bas_loc, mol_loc, true);
            val comp_loc = new ComputePlaceNewFuture(den_loc, twoE_loc);
            at(computeInst) { computeInst(PIdx(0)++) = comp_loc; };
            
            Console.OUT.println("Starting computation..." + here);

            val F1 = future(G) { 
                       var myG:Int; 
                       atomic myG = G(0)++;
                       return myG;
            };

            myG = F1.force();

            for(a=0; a<noOfAtoms; a++) {
              val aFunc = mol_loc.getAtom(a).getBasisFunctions();
              naFunc = aFunc.size();
              // basis functions on a
              for(i=0; i<naFunc; i++) {
                val iaFunc = aFunc.get(i);

                // center b
                for(b=0; b<=a; b++) {
                    val bFunc = mol_loc.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            val cFunc = mol_loc.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                val kcFunc = cFunc.get(k);

                                // center d
                                for(d=0; d<=c; d++) {
                                    val dFunc = mol_loc.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        val ldFunc = dFunc.get(l);

                                        if (L == myG) {
                                          val F2 = future(G) {
                                                 var myG:Int;
                                                 atomic myG = G(0)++;
                                                 return myG;
                                          };

            	                          comp_loc.compute2EAndRecord(iaFunc, jbFunc, kcFunc, ldFunc); 

                                          myG = F2.force();
                                        } // end if

                                        L++;
				    } // end l
				} // center d
                            } // end k
                        } // center c
                    } // end j
                } // center b
              } // end i
            } // center a
          } // at
        } // foreach
        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds");
    
        timer.start(1);

        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  /**
                  val x_loc = x; val y_loc = y;
                  val jVal = at(comp_loc) { comp_loc.getJMat().getMatrix()(x_loc, y_loc) };
                  val kVal = at(comp_loc) { comp_loc.getKMat().getMatrix()(x_loc, y_loc) };

                  gMatrix(x,y) += jVal - (0.25*kVal);
                  **/

                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for

        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds");
    }


    /** find unique elements and mark the onces that are not */
    /** 8 => is the level of integral symmetry, given (i,j|k.l)
        there are 8 combinations that are unique */
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

    /** Compute class for the old code */
    class ComputeOld {
        var computing:Boolean = false;

        var i:Int, j:Int, k:Int, l:Int;

        val twoEI:TwoElectronIntegrals!;
        val jMat:Matrix!, kMat:Matrix!;
        val density:Density!;
        val idx_loc:Rail[Int]!;
        val jdx_loc:Rail[Int]!;
        val kdx_loc:Rail[Int]!;
        val ldx_loc:Rail[Int]!;
        val validIdx:Rail[Boolean]!;

        val jMatrix:Array[Double]{rank==2};
        val kMatrix:Array[Double]{rank==2};
        val dMatrix:Array[Double]{rank==2};
    
        public def this(te:TwoElectronIntegrals, den:Density) {
            twoEI = te as TwoElectronIntegrals!;
            density = den as Density!;

            idx_loc = Rail.make[Int](8);
            jdx_loc = Rail.make[Int](8);
            kdx_loc = Rail.make[Int](8);
            ldx_loc = Rail.make[Int](8);
 
            validIdx = Rail.make[Boolean](8, (Int)=>true);

            val N = density.getRowCount();

            jMat = new Matrix(N) as Matrix!;
            kMat = new Matrix(N) as Matrix!;

            jMatrix = jMat.getMatrix();
            kMatrix = kMat.getMatrix();
            dMatrix = density.getMatrix();

            jMat.makeZero();
            kMat.makeZero();
        }

        public def compute() {
            val twoEIntVal = twoEI.compute2E(i, j, k, l);

            for(var m:Int=0; m<8; m++) validIdx(m) = true;

            setJKMatrixElements(i, j, k, l, twoEIntVal);

            // if not special case
            if ((i|j|k|l) != 0) {
               // else this is symmetry unique integral, so need to
               // use this value for all 8 combinations
               // (if unique)
               ldx_loc(0) = l; ldx_loc(1) = l; kdx_loc(2) = l; kdx_loc(3) = l;
               idx_loc(4) = l; idx_loc(5) = l; jdx_loc(6) = l; jdx_loc(7) = l;

               // filter unique elements
               filterUniqueElements();

               // and evaluate them
               for(var m:Int=1; m<8; m++) {
                   if (validIdx(m)) {
                       setJKMatrixElements(idx_loc(m), jdx_loc(m), kdx_loc(m), ldx_loc(m), twoEIntVal);
                   } // end if
               } // end m
            } // end if

            atomic computing = false;
        }

        /** find unique elements and mark the onces that are not */
        private def filterUniqueElements() : void {
           var i:Int, j:Int, k:Int, l:Int, m:Int, n:Int;

           for(m=0; m<8; m++) {
              i = idx_loc(m); j = jdx_loc(m); k = kdx_loc(m); l = ldx_loc(m);
              for(n=m+1; n<8; n++) {
                  if (i==idx_loc(n) && j==jdx_loc(n) && k==kdx_loc(n) && l==ldx_loc(n))
                      validIdx(n) = false;
              } // end for
           } // end for
        }

        /** Set the J and K value for a given combination */
        private def setJKMatrixElements(i:Int, j:Int, k:Int, l:Int, twoEIntVal:Double) : void {
           val v1 = dMatrix(k,l) * twoEIntVal;
           val v2 = dMatrix(i,j) * twoEIntVal;
           val v3 = dMatrix(j,l) * twoEIntVal;
           val v4 = dMatrix(j,k) * twoEIntVal;
           val v5 = dMatrix(i,l) * twoEIntVal;
           val v6 = dMatrix(i,k) * twoEIntVal;

           jMatrix(i,j) += v1;
           jMatrix(k,l) += v2;
           kMatrix(i,k) += v3;
           kMatrix(i,l) += v4;
           kMatrix(j,k) += v5;
           kMatrix(j,l) += v6;
        }

        public def setValue(i:Int, j:Int, k:Int, l:Int, 
                            idx:Rail[Int]!, jdx:Rail[Int]!, kdx:Rail[Int]!, ldx:Rail[Int]!) : Boolean {
            if (computing) return false;

            this.i = i;
            this.j = j;
            this.k = k;
            this.l = l;

            for(var m:Int=0; m<8; m++) { 
               idx_loc(m) = idx(m);
               jdx_loc(m) = jdx(m);
               kdx_loc(m) = kdx(m);
               ldx_loc(m) = ldx(m);
            } // end for

            atomic computing = true;

            return true;
        }

        public def getJMat() = jMat;
        public def getKMat() = kMat;
    }


    /** Compute class for the new code */
    class ComputeNew {
        var computing:Boolean = false;

        var i:ContractedGaussian!, j:ContractedGaussian!, 
            k:ContractedGaussian!, l:ContractedGaussian!;

        val twoEI:TwoElectronIntegrals!;
        val shellList:ShellList;
        val jMat:Matrix!, kMat:Matrix!;
        val density:Density;

        public def this(te:TwoElectronIntegrals, sh:ShellList, den:Density) { 
            twoEI = te as TwoElectronIntegrals!;
            shellList = sh;
            density = den;

            val N = density.getRowCount();

            jMat = new Matrix(N) as Matrix!;
            kMat = new Matrix(N) as Matrix!;

            jMat.makeZero();
            kMat.makeZero();
        }

        public def compute() {
            twoEI.compute2EAndRecord(i as ContractedGaussian!, j as ContractedGaussian!, 
                                     k as ContractedGaussian!, l as ContractedGaussian!,  
                                     shellList as ShellList!, 
                                     jMat as Matrix!, kMat as Matrix!, density as Density!);
            atomic computing = false;
        }
        
        public def setValue(i:ContractedGaussian!, j:ContractedGaussian!,
                            k:ContractedGaussian!, l:ContractedGaussian!) : Boolean {
            if (computing) return false;

            this.i = i;
            this.j = j;
            this.k = k;
            this.l = l;
            atomic computing = true;

            return true;
        }

        public def getJMat() = jMat;
        public def getKMat() = kMat;
    }

    /** Compute class for the old code - multi place version */
    class ComputePlaceOld {
      
        val thePlace:Place;

        val computeInst = Rail.make[ComputeOld!](Runtime.INIT_THREADS);

        val idx = Rail.make[Int](8);
        val jdx = Rail.make[Int](8);
        val kdx = Rail.make[Int](8);
        val ldx = Rail.make[Int](8);

        val density:Density!;

        public def this(mol:Molecule[QMAtom], bfs:BasisFunctions, den:Density, tp:Place) {
            thePlace = tp;

            val mol_loc = new Molecule[QMAtom]() as Molecule[QMAtom]!;
            val nAtoms  = at(mol) { mol.getNumberOfAtoms() };
 
            for(var i:Int=0; i<nAtoms; i++) {
                val i_loc = i;
                val sym  = at(mol) { mol.getAtom(i_loc).symbol };
                val x    = at(mol) { mol.getAtom(i_loc).centre.i };
                val y    = at(mol) { mol.getAtom(i_loc).centre.j };
                val z    = at(mol) { mol.getAtom(i_loc).centre.k };
   
                mol_loc.addAtom(new QMAtom(sym, new Point3d(x, y, z)));
            } // end for

            // Console.OUT.println("\tStart making local Molecule and BasisFunctions @ " + here);

            val basisName = at(bfs) { bfs.getBasisName() };
            val bas_loc = new BasisFunctions(mol_loc, basisName, "basis");

            // Console.OUT.println("\tDone making local Molecule and BasisFunctions @ " + here);

            // Console.OUT.println("\tMake local copy of density @ " + here);

            val N = at(den) { den.getRowCount() };
            val den_loc = new Density(N);
            density = den_loc;
            val den_loc_mat = den_loc.getMatrix();
            val den_rem = at(den) { den.getValRail() };
            var i:Int, j:Int, ii:Int;

            ii = 0 ;  
            for(i=0; i<N; i++) 
              for(j=0; j<N; j++) 
                  den_loc_mat(i, j) = den_rem(ii++);

            // Console.OUT.println("\tDone making local copy of density @ " + here);

            for(i=0; i<Runtime.INIT_THREADS; i++) {
                computeInst(i) = new ComputeOld(new TwoElectronIntegrals(bas_loc, mol_loc, true), den_loc);
            } // end for

            // Console.OUT.println("\tDone initing tasks @ " + here);
        }

        public def setValue(i:Int, j:Int, k:Int, l:Int) : Boolean {
            idx(0) = i; jdx(1) = i; jdx(2) = i; idx(3) = i;
            kdx(4) = i; ldx(5) = i; kdx(6) = i; ldx(7) = i;

            jdx(0) = j; idx(1) = j; idx(2) = j; jdx(3) = j;
            ldx(4) = j; kdx(5) = j; ldx(6) = j; kdx(7) = j;

            kdx(0) = k; kdx(1) = k; ldx(2) = k; ldx(3) = k;
            jdx(4) = k; jdx(5) = k; idx(6) = k; idx(7) = k;

            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
               val setIt = computeInst(ix).setValue(i, j, k, l, idx, jdx, kdx, ldx);

               if (setIt) {
                  val ix_loc = ix;
                  async computeInst(ix_loc).compute();
                  return true;
               } // end if
            } // end for

            return false;
        }

        public def getJMat() : Matrix! {
            val N = density.getRowCount();
            var jM:Matrix! = new Matrix(N) as Matrix!;

            jM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               jM = jM.add(computeInst(i).getJMat());
            } //  end for

            return jM;             
        }

        public def getKMat() : Matrix! {
            val N = density.getRowCount();
            var kM:Matrix! = new Matrix(N) as Matrix!;

            kM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               kM = kM.add(computeInst(i).getKMat());
            } //  end for

            return kM;
        }

        // TODO following two methods should not be necessary once XTENLANG-787 is resolved
        public def getJMatVal() : ValRail[Double]! {
           return getJMat().getValRail();
        }

        public def getKMatVal() : ValRail[Double]! {
           return getKMat().getValRail();
        }
    }

    /** Compute class for the old code - multi place version , using a task pool */
    class ComputePlaceOldPool extends ComputePlaceOld {
        val taskPool:Stack[BlockIndices]!;
        var noMoreTasks:Boolean;

        public def this(mol:Molecule[QMAtom], bfs:BasisFunctions, den:Density, tp:Place) {
           super(mol, bfs, den, tp);

           taskPool = new Stack[BlockIndices]() as Stack[BlockIndices]!;
           noMoreTasks = false;
        }

        public def pushBlockIdx(bIdx:BlockIndices) {
           taskPool.push(bIdx);
        }

        public def start() {
           outer: while(true) {
               for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                   if (noMoreTasks && (taskPool.size()==0)) break outer;

                   await (computeInst(ix).computing == false);

                   if (taskPool.size() == 0) continue;

                   val task = taskPool.pop();

                   idx(0) = task.i; jdx(1) = task.i; jdx(2) = task.i; idx(3) = task.i;
                   kdx(4) = task.i; ldx(5) = task.i; kdx(6) = task.i; ldx(7) = task.i;

                   jdx(0) = task.j; idx(1) = task.j; idx(2) = task.j; jdx(3) = task.j;
                   ldx(4) = task.j; kdx(5) = task.j; ldx(6) = task.j; kdx(7) = task.j;

                   kdx(0) = task.k; kdx(1) = task.k; ldx(2) = task.k; ldx(3) = task.k;
                   jdx(4) = task.k; jdx(5) = task.k; idx(6) = task.k; idx(7) = task.k;

                   computeInst(ix).setValue(task.i, task.j, task.k, task.l, idx, jdx, kdx, ldx);
                   val ix_loc = ix;
                   async computeInst(ix_loc).compute();
               } // end for
           } // end while
        }
    }

    /** Compute class for the old code - multi place version , just a direct compute wrapper */
    class ComputePlaceOldDirect extends ComputePlaceOld {
        public def this(mol:Molecule[QMAtom], bfs:BasisFunctions, den:Density, tp:Place) {
              super(mol, bfs, den, tp);
        }

        public def compute(i:Int, j:Int, k:Int, l:Int) {
              // TODO: this is actually handling only one thread
              for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                   await (computeInst(ix).computing == false);

                   computeInst(ix).setValue(i, j, k, l, idx, jdx, kdx, ldx);
                   val ix_loc = ix;
                   async computeInst(ix_loc).compute();
                   if (ix == 0) break;
              } // end for
        }

        public def compute(start:Int, end:Int) {
              val N = density.getRowCount();
              var i:Int, j:Int, k:Int, l:Int, m:Int, ij:Int, kl:Int;

              for(i=start; i<end; i++) {
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
                              compute(i, j, k, l);
                          }
                       }
                   }
                }
              }
        }
    }

    /** Compute class for the new code - multi place version , just a direct compute wrapper */
    class ComputePlaceNewDirect extends ComputePlaceNew {

        val molecule:Molecule[QMAtom]!;

        public def this(mol:Molecule[QMAtom], bfs:BasisFunctions, den:Density, tp:Place) {
            super(mol, bfs, den, tp);

            molecule = mol_loc;
        }

        public def compute(i:ContractedGaussian, j:ContractedGaussian, 
                           k:ContractedGaussian, l:ContractedGaussian) {
            // TODO: this is actually handling only one thread
            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                await (computeInst(ix).computing == false);

                computeInst(ix).setValue(i as ContractedGaussian!, j as ContractedGaussian!, 
                                         k as ContractedGaussian!, l as ContractedGaussian!);
                val ix_loc = ix;
                async computeInst(ix_loc).compute();
                if (ix == 0) break;
            } // end for
        }

        public def compute(start:Int, end:Int) {
            var i:Int, j:Int, k:Int, l:Int;
            var a:Int, b:Int, c:Int, d:Int;
            var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
                kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};
            var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int;
            var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
                bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
                cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
                dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};

            val noOfAtoms = molecule.getNumberOfAtoms();

            for(a=start; a<end; a++) {
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

                                       // TODO: 
                                       compute(iaFunc, jbFunc, kcFunc, ldFunc);
                                   } // end l
                               } // center d
                           } // end k
                       } // center c
                   } // end j
               } // center b
             } // end i
          } // center a
        }
    } 

    /** Compute class for the new code - multi place version */
    class ComputePlaceNew {
        val thePlace:Place;

        val computeInst = Rail.make[ComputeNew!](Runtime.INIT_THREADS);

        val density:Density!;

        val bas_loc:BasisFunctions!;
        val mol_loc:Molecule[QMAtom]!;

        public def this(mol:Molecule[QMAtom], bfs:BasisFunctions, den:Density, tp:Place) {
            thePlace = tp;

            val mol_loc = new Molecule[QMAtom]() as Molecule[QMAtom]!;
            val nAtoms  = at(mol) { mol.getNumberOfAtoms() };
 
            for(var i:Int=0; i<nAtoms; i++) {
                val i_loc = i;
                val sym  = at(mol) { mol.getAtom(i_loc).symbol };
                val x    = at(mol) { mol.getAtom(i_loc).centre.i };
                val y    = at(mol) { mol.getAtom(i_loc).centre.j };
                val z    = at(mol) { mol.getAtom(i_loc).centre.k };
   
                mol_loc.addAtom(new QMAtom(sym, new Point3d(x, y, z)));
            } // end for

            this.mol_loc = mol_loc;

            // Console.OUT.println("\tStart making local Molecule and BasisFunctions @ " + here);

            val basisName = at(bfs) { bfs.getBasisName() };
            val bas_loc = new BasisFunctions(mol_loc, basisName, "basis");

            this.bas_loc = bas_loc;

            // Console.OUT.println("\tDone making local Molecule and BasisFunctions @ " + here);

            // Console.OUT.println("\tMake local copy of density @ " + here);

            val N = at(den) { den.getRowCount() };
            val den_loc = new Density(N);
            density = den_loc;
            val den_loc_mat = den_loc.getMatrix();
            val den_rem = at(den) { den.getValRail() };
            var i:Int, j:Int, ii:Int;

            ii = 0 ;
            for(i=0; i<N; i++)
              for(j=0; j<N; j++)
                  den_loc_mat(i, j) = den_rem(ii++);

            // Console.OUT.println("\tDone making local copy of density @ " + here);

            for(i=0; i<Runtime.INIT_THREADS; i++) {
                computeInst(i) = new ComputeNew(new TwoElectronIntegrals(bas_loc, mol_loc, true), bas_loc.getShellList(), den_loc);          
            } // end for

            // Console.OUT.println("\tDone initing tasks @ " + here);
        }

        public def setValue(a:Int, b:Int, c:Int, d:Int, i:Int, j:Int, k:Int, l:Int) : Boolean {
            val iFunc = mol_loc.getAtom(a).getBasisFunctions().get(i);
            val jFunc = mol_loc.getAtom(b).getBasisFunctions().get(j);
            val kFunc = mol_loc.getAtom(c).getBasisFunctions().get(k);
            val lFunc = mol_loc.getAtom(d).getBasisFunctions().get(l);

            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
               val setIt = computeInst(ix).setValue(iFunc, jFunc, kFunc, lFunc);

               if (setIt) {
                  val ix_loc = ix;
                  async computeInst(ix_loc).compute();
                  return true;
               } // end if
            } // end for

            return false;
        }

        public def getJMat() : Matrix! {
            val N = density.getRowCount();
            var jM:Matrix! = new Matrix(N) as Matrix!;

            jM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               jM = jM.add(computeInst(i).getJMat());
            } //  end for

            return jM;             
        }

        public def getKMat() : Matrix! {
            val N = density.getRowCount();
            var kM:Matrix! = new Matrix(N) as Matrix!;

            kM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               kM = kM.add(computeInst(i).getKMat());
            } //  end for

            return kM;
        }

        // TODO following two methods should not be necessary once XTENLANG-787 is resolved
        public def getJMatVal() : ValRail[Double]! {
           return getJMat().getValRail();
        }

        public def getKMatVal() : ValRail[Double]! {
           return getKMat().getValRail();
        }
    }

    class ComputePlaceNewFuture {
         val density:Density!;
         val twoEI:TwoElectronIntegrals!;

         val jM:Matrix!;
         val kM:Matrix!;
         val shellList:ShellList!; 

         public def this(den:Density!, twoEI:TwoElectronIntegrals!) { 
             this.density = den;
             this.twoEI   = twoEI;

             shellList = twoEI.getBasisFunctions().getShellList();

             val N = den.getRowCount();

             jM = new Matrix(N) as Matrix!;
             kM = new Matrix(N) as Matrix!;
         }

         public def compute2EAndRecord(iaFunc:ContractedGaussian!, jbFunc:ContractedGaussian!, 
                                       kcFunc:ContractedGaussian!, ldFunc:ContractedGaussian!) {
             twoEI.compute2EAndRecord(iaFunc, jbFunc, kcFunc, ldFunc, shellList, jM, kM, density);
         }

         public def getJMat() : Matrix! {
             return jM;
         }
         
         public def getKMat() : Matrix! {
             return kM;
         }

         // TODO following two methods should not be necessary once XTENLANG-787 is resolved
         public def getJMatVal() : ValRail[Double]! {
            return getJMat().getValRail();
         }

         public def getKMatVal() : ValRail[Double]! {
            return getKMat().getValRail();
         }
    }

    static struct BlockIndices(i:Int, j:Int, k:Int, l:Int) {
         public def this(i:Int, j:Int, k:Int, l:Int) {
             property(i, j, k, l);
         }
    }
}

