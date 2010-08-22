/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import x10.util.*;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;

/**
 * GMatrix.x10
 *
 * GMatrix in HF calculation
 *
 * @author: V.Ganesh
 */
public class GMatrix extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    /** top level method to form the G Matrix, depending on gMatType appropriate functions are called */
    public def compute(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, density:Density!, gMatType:Int) : void {
       val timer = new Timer(1);

       timer.start(0);
       switch(gMatType) {
       case 0:
           Console.OUT.println("   GMatrix.computeDirectSerialNew: " + gMatType);
           computeDirectSerialNew(bfs, mol, density); 
           break;
       case 1:
           Console.OUT.println("   GMatrix.computeDirectLowMemNewNoAtomic: " + gMatType);
           computeDirectLowMemNewNoAtomic(bfs, mol, density); 
           break;
       case 2:
           Console.OUT.println("   GMatrix.computeDirectNewMultiPlaceNoAtomic: " + gMatType);
           computeDirectNewMultiPlaceNoAtomic(bfs, mol, density);
           break;
       case 3:
           Console.OUT.println("   GMatrix.computeDirectNewMultiPlaceStatic: " + gMatType);
           computeDirectNewMultiPlaceStatic(bfs, mol, density);
           break;
       case 4:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceNewFuture: " + gMatType);
           computeDirectMultiPlaceNewFuture(bfs, mol, density);
           break;
       case 5:
       default:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceShellLoop: " + gMatType);
           computeDirectMultiPlaceShellLoop(bfs, mol, density);
           break;
       } // end switch .. case
       timer.stop(0);
       Console.OUT.println ("    Time to construct GMatrix: " + (timer.total(0) as Double) / 1e9 + " seconds");
    }

    private def computeDirectSerialNew(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val jMat = new Matrix(N);
        val kMat = new Matrix(N);

        jMat.makeZero();
        kMat.makeZero();

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();

        var i:Int, j:Int, k:Int, l:Int;

        val shellList = bfs.getShellList();
        val noOfBasisFunctions = bfs.getBasisFunctions().size();
        

        val noOfAtoms = mol.getNumberOfAtoms();
        var a:Int, b:Int, c:Int, d:Int;
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)},
            dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)},
            kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};

        val maxam = shellList.getMaximumAngularMomentum();
        val twoE = new TwoElectronIntegrals(maxam);

        // center a
        for(a=0; a<noOfAtoms; a++) {
            aFunc = mol.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
                iaFunc = aFunc.get(i);

                // center b
                for(b=0; b<=a; b++) {
                    bFunc = mol.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        jbFunc = bFunc.get(j);

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            cFunc = mol.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                kcFunc = cFunc.get(k);

                                // center d
                                for(d=0; d<=c; d++) {
                                    dFunc = mol.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        ldFunc = dFunc.get(l);

                                        // Console.OUT.println(a + ", " + b + ", " + c + ", " + d);

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
        val gMatrix = getMatrix();
        finish foreach((x,y) in gMatrix.region)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    private def computeDirectLowMemNewNoAtomic(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, 
                                               density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        val shellList = bfs.getShellList();
        val noOfBasisFunctions = bfs.getBasisFunctions().size();

        val noOfAtoms = mol.getNumberOfAtoms();

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

        val maxam = shellList.getMaximumAngularMomentum();
        for(i=0; i<Runtime.INIT_THREADS; i++) {
            computeInst(i) = new ComputeNew(new TwoElectronIntegrals(maxam),
                                            shellList, density);
        } // end for

        // center a
        finish {
          for(a=0; a<noOfAtoms; a++) {
            aFunc = mol.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
               iaFunc = aFunc.get(i);

               // center b
               for(b=0; b<=a; b++) {
                   bFunc = mol.getAtom(b).getBasisFunctions();
                   nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(j=0; j<nbFunc; j++) {
                       jbFunc = bFunc.get(j);

                       // center c
                       for(c=0; c<noOfAtoms; c++) {
                           cFunc = mol.getAtom(c).getBasisFunctions();
                           ncFunc = cFunc.size();
                           // basis functions on c
                           for(k=0; k<ncFunc; k++) {
                               kcFunc = cFunc.get(k);

                               // center d
                               for(d=0; d<=c; d++) {
                                   dFunc = mol.getAtom(d).getBasisFunctions();
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
        val gMatrix = getMatrix();
        for(var idx:Int=0; idx<Runtime.INIT_THREADS; idx++) { 
             val jMatrix = computeInst(idx).getJMat().getMatrix();
             val kMatrix = computeInst(idx).getKMat().getMatrix();
             
             finish foreach((x,y) in gMatrix) {
                   gMatrix(x,y) += jMatrix(x,y) - (0.25*kMatrix(x,y));     
             } // finish
        } // end for
    }

    private def computeDirectNewMultiPlaceNoAtomic(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, 
                                                   density:Density!) : void {
        makeZero();

        val noOfAtoms = mol.getNumberOfAtoms();

        val timer = new Timer(3);

        timer.start(0);
        val basisName = bfs.getBasisName();
        val computeInst = DistArray.make[ComputePlaceNew](Dist.makeUnique(Place.places), ((p) : Point) => new ComputePlaceNewDirect(mol, basisName, density));
        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        finish {
            // center a
            for ((a) in 0..(noOfAtoms-1)) {
                val aFunc = mol.getAtom(a).getBasisFunctions();
                val naFunc = aFunc.size();
                // basis functions on a
                for ((i) in 0..(naFunc-1)) {
                    val iaFunc = aFunc.get(i);

                    // center b
                    for ((b) in 0..a) {
                        val bFunc = mol.getAtom(b).getBasisFunctions();
                        val nbFunc = (b<a) ? bFunc.size() : i+1;
                        // basis functions on b
                        for ((j) in 0..(nbFunc-1)) {
                            val jbFunc = bFunc.get(j);

                            // center c
                            for ((c) in 0..(noOfAtoms-1)) {
                                val cFunc = mol.getAtom(c).getBasisFunctions();
                                val ncFunc = cFunc.size();
                                // basis functions on c
                                for ((k) in 0..(ncFunc-1)) {
                                    val kcFunc = cFunc.get(k);

                                    // center d
                                    for ((d) in 0..c) {
                                       val dFunc = mol.getAtom(d).getBasisFunctions();
                                       val ndFunc = (d<c) ? dFunc.size() : k+1;
                                        // basis functions on d
                                        for ((l) in 0..(ndFunc-1)) {
                                            val ldFunc = dFunc.get(l);
                                            var setIt:Boolean = false;

                                            outer: while(!setIt) {
                                                for(p in computeInst) {
                                                    setIt = at(computeInst.dist(p)) { computeInst(p).setValue(a, b, c, d, i, j, k, l) };

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
        val N = density.getRowCount();
        val gMatrix = getMatrix();
        for(p in computeInst) {
             val jVal = at(computeInst.dist(p)) { computeInst(p).getJMatVal() };
             val kVal = at(computeInst.dist(p)) { computeInst(p).getKMatVal() };

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

    private def computeDirectNewMultiPlaceStatic(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, 
                                                 density:Density!) : void {
        val N = density.getRowCount();

        makeZero();

        var i:Int, j:Int;

        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceNewDirect](nPlaces);

        val timer = new Timer(3);

        timer.start(0);
        i = 0;
        val basisName = bfs.getBasisName();
        for(place in Place.places) {
           computeInst(i) = at(place) { return new ComputePlaceNewDirect(mol, basisName, density); };
           i++;
        }

        val workPerPlace = Rail.make[Int](nPlaces, (Int)=>0);
        
        i = mol.getNumberOfAtoms();
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
        val gMatrix = getMatrix();
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
    private def computeDirectMultiPlaceNewFuture(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, density:Density!) : void {
        // init counters
        G(0) = 0; PIdx(0) = 0;

        // init other variables
        val N = density.getRowCount();

        makeZero();

        val shellList = bfs.getShellList();

        val noOfAtoms = mol.getNumberOfAtoms();
        val nPlaces = Place.places.length;
        val computeInst = Rail.make[ComputePlaceNewFuture](nPlaces);
  
        val timer = new Timer(2);

        timer.start(0);

        val basisName = bfs.getBasisName();
    
        // center a
        finish foreach(place in Place.places) {
          async at(place) {
            var myG:Int = 0;
            var L:Int = 0;

            var i:Int, j:Int, k:Int, l:Int;
            var a:Int, b:Int, c:Int, d:Int;
            var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;

            // TODO: better way to pass global data instead?
            // make local copies of data 
            val mol_loc = new Molecule[QMAtom]();

            for((i) in 0..(nAtoms-1)) {
                val sym  = at(mol) { mol.getAtom(i).symbol };
                val centre = at(mol) {mol.getAtom(i).centre };
                mol_loc.addAtom(new QMAtom(sym, centre));
            }

            val den_loc = new Density(density);

            val bfs = new BasisFunctions(mol_loc, basisName, "basis");
            val comp_loc = new ComputePlaceNewFuture(bfs, den_loc);
            at(computeInst) { computeInst(PIdx(0)++) = comp_loc; };

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
        val gMatrix = getMatrix();
        for(comp_loc in computeInst) {
             val jVal = at(comp_loc) { comp_loc.getJMatVal() };
             val kVal = at(comp_loc) { comp_loc.getKMatVal() };

             var ii:Int=0;
             for(var x:Int=0; x<N; x++) {
               for(var y:Int=0; y<N; y++) {
                  gMatrix(x,y) += jVal(ii) - (0.25*kVal(ii));
                  ii++;
               } // end for
             } // end for
        } // end for

        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds");
    }

    private def computeDirectMultiPlaceShellLoop(bfs:BasisFunctions!, mol:Molecule[QMAtom]!, density:Density!) : void {
        makeZero();

        val shellList = bfs.getShellList();
        val shellPairs = shellList.getShellPairs();

        val timer = new Timer(3);

        timer.start(0);
        val nPairs = shellPairs.length();
        val basisName = bfs.getBasisName();
        val computeInst = DistArray.make[ComputePlaceNewDirect](Dist.makeUnique(Place.places), ((p) : Point) => new ComputePlaceNewDirect(mol, basisName, density));

        timer.stop(0);
        Console.OUT.println("\tTime for setting up place(s) with initial data: " + (timer.total(0) as Double) / 1e9 + " seconds");

        val nPlaces = Place.places.length;
        val workPerPlace = nPairs / nPlaces;
        val remainder = nPairs % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);
        Console.OUT.println("\tWork units per place: " + workPerPlace + " remainder " + remainder);

        timer.start(1);

        finish ateach ((placeId) in computeInst) {
            val start = placeId < remainder ? (placeId * (workPerPlace + 1)) : (firstChunk + (placeId-remainder) * workPerPlace);
            val end = start + workPerPlace + (placeId < remainder ? 1 : 0);
            computeInst(placeId).computeShells(start, end, nPairs);
        }

        timer.stop(1);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(1) as Double) / 1e9 + " seconds");

        timer.start(2);
        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        val gMatrix = getMatrix();
        val N = density.getRowCount();
        for(p in computeInst) {
            val gVal = at(computeInst.dist(p)) { computeInst(p).getGMatVal() };

            // add place contribution to gMatrix
            var ii:Int = 0;
            for(var x:Int=0; x<N; x++) {
                for(var y:Int=0; y<N; y++) {
                  gMatrix(x,y) += gVal(ii);
                  ii++;
               } // end for
             } // end for

        } // end for
        timer.stop(2);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(2) as Double) / 1e9 + " seconds");

    }

    /** Compute class for the new code */
    class ComputeNew {
        var computing:Boolean = false;

        var i:ContractedGaussian!, j:ContractedGaussian!, 
            k:ContractedGaussian!, l:ContractedGaussian!;

        val twoEI:TwoElectronIntegrals!;
        val shellList:ShellList{self.at(this)};
        val jMat:Matrix!, kMat:Matrix!;
        val density:Density!;

        public def this(te:TwoElectronIntegrals!, sh:ShellList, den:Density!) { 
            twoEI = te;
            shellList = sh;
            density = den;

            val N = density.getRowCount();

            jMat = new Matrix(N);
            kMat = new Matrix(N);

            jMat.makeZero();
            kMat.makeZero();
        }

        public def compute() {
            twoEI.compute2EAndRecord(i, j, 
                                     k, l,  
                                     shellList, 
                                     jMat, kMat, density);
            atomic computing = false;
        }

        public def computeSingle(i:ContractedGaussian!, j:ContractedGaussian!,
                                 k:ContractedGaussian!, l:ContractedGaussian!) {
            twoEI.compute2EAndRecord(i, j, k, l, 
                                     shellList,
                                     jMat, kMat, density);
        }

        public def computeSingle2(i:ContractedGaussian!, j:ContractedGaussian!,
                                  k:ContractedGaussian!, l:ContractedGaussian!, 
                                  radiusABSquared:Double,
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int, abLim:Int) {
            twoEI.compute2EAndRecord2(i, j, k, l,
                                      shellList,
                                      jMat, kMat, density,
                                      radiusABSquared,
                                      aAng, bAng, cAng, dAng, angMomAB,
                                      aStrt, bStrt, cStrt, dStrt,
                                      aLim, bLim, abLim);
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

    /** Compute class for the new code - multi place version , just a direct compute wrapper */
    class ComputePlaceNewDirect extends ComputePlaceNew {
        val gMatrixContribution : Matrix!;

        public def this(mol:Molecule[QMAtom], basisName:String, den:Density) {
            super(mol, basisName, den);

            val N = den.getRowCount();
            gMatrixContribution = new Matrix(N);
        }

        public def compute(i:ContractedGaussian!, j:ContractedGaussian!, 
                           k:ContractedGaussian!, l:ContractedGaussian!) {
            // TODO: this is actually handling only one thread
            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                await (computeInst(ix).computing == false);

                computeInst(ix).setValue(i, j, k, l);
                val ix_loc = ix;
                async computeInst(ix_loc).compute();
                if (ix == 0) break;
            } // end for
        }

        public def compute(start:Int, end:Int) {
            val noOfAtoms = mol_loc.getNumberOfAtoms();

            for(var a:Int=start; a<end; a++) {
             val aFunc = mol_loc.getAtom(a).getBasisFunctions();
             val naFunc = aFunc.size();
             // basis functions on a
             for(var i:Int=0; i<naFunc; i++) {
               val iaFunc = aFunc.get(i);

               // center b
               for(var b:Int=0; b<=a; b++) {
                   val bFunc = mol_loc.getAtom(b).getBasisFunctions();
                   val nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(var j:Int=0; j<nbFunc; j++) {
                       val jbFunc = bFunc.get(j);

                       // center c
                       for(var c:Int=0; c<noOfAtoms; c++) {
                           val cFunc = mol_loc.getAtom(c).getBasisFunctions();
                           val ncFunc = cFunc.size();
                           // basis functions on c
                           for(var k:Int=0; k<ncFunc; k++) {
                               val kcFunc = cFunc.get(k);

                               // center d
                               for(var d:Int=0; d<=c; d++) {
                                   val dFunc = mol_loc.getAtom(d).getBasisFunctions();
                                   val ndFunc = (d<c) ? dFunc.size() : k+1;
                                   // basis functions on d
                                   for(var l:Int=0; l<ndFunc; l++) {
                                       val ldFunc = dFunc.get(l);

                                       // TODO: 
                                       // Console.OUT.println(a + ", " + b + ", " + c + ", " + d + " | " + i + ", " + j + ", " + k + ", " + l);
                                       computeInst(0).computeSingle(iaFunc, jbFunc, kcFunc, ldFunc);
                                   } // end l
                               } // center d
                           } // end k
                       } // center c
                   } // end j
               } // center b
             } // end i
          } // center a
        }

        public def computeShells(startShell:Int, endShell:Int, nPairs:Int) {
               val bfs = computeInst(0).shellList.getShellPrimitives();
               val shellPairs = computeInst(0).shellList.getShellPairs();

               for(var i:Int=startShell; i<endShell; i++) {
                  val a = shellPairs(i).first;
                  val b = shellPairs(i).second;

                  val aFunc = bfs(a) as ContractedGaussian!;
                  val bFunc = bfs(b) as ContractedGaussian!;

                  val aStrt = aFunc.getIntIndex();
                  val bStrt = bFunc.getIntIndex();
                  val aAng  = aFunc.getMaximumAngularMomentum();
                  val bAng  = bFunc.getMaximumAngularMomentum();

                  val aa = aStrt + aAng;
                  val bb = bStrt + bAng;

                  if (aa < bb) continue;

                  val angMomAB = aAng + bAng;
                  val aLim = ((aAng+1)*(aAng+2)/2);
                  val bLim = ((bAng+1)*(bAng+2)/2);
                  val abLim = aLim * bLim;

                  val radiusABSquared = aFunc.distanceSquaredFrom(bFunc);

                  for(var j:Int=0; j<nPairs; j++) {
                     val c = shellPairs(j).first;
                     val cFunc = bfs(c) as ContractedGaussian!;
                     val cStrt = cFunc.getIntIndex();
                     val cAng  = cFunc.getMaximumAngularMomentum();
                     val cc = cStrt + cAng;

                     if (aa < cc) continue;

                     val d = shellPairs(j).second;
                     val dFunc = bfs(d) as ContractedGaussian!;
                     val dStrt = dFunc.getIntIndex();
                     val dAng  = dFunc.getMaximumAngularMomentum();
                     val dd = dStrt + dAng;

                     if (cc < dd) continue;

                     // if (bb < cc && bb < dd) continue;
                     // val skp = (bb < cc && bb < dd && aa == cc);

                     // Console.OUT.println(aa + ", " + bb + ", " + cc + ", " + dd + " : (" + a + "," + b + "," + c + "," + d + ") " + skp);
                     // twoE.compute2EAndRecord(aFunc, bFunc, cFunc, dFunc, shellList, jMat, kMat, density);
                     computeInst(0).computeSingle2(aFunc, bFunc, cFunc, dFunc,
                                    radiusABSquared, aAng, bAng, cAng, dAng, angMomAB,
                                    aStrt, bStrt, cStrt, dStrt, aLim, bLim, abLim);
                  } // end for
               } // end for

            // calculate this place's contribution to G Matrix
            val jMat = getJMat().getMatrix();
            val kMat = getKMat().getMatrix();
            val gMat = gMatrixContribution.getMatrix();
            for (p in jMat.region) {
                gMat(p) = jMat(p) - 0.25 * kMat(p);
            }
        }

        // TODO following method should not be necessary once XTENLANG-787 is resolved
        public def getGMatVal() : ValRail[Double]! {
           return gMatrixContribution.getValRail();
        }
    } 

    /** Compute class for the new code - multi place version */
    class ComputePlaceNew {
        val computeInst = Rail.make[ComputeNew!](Runtime.INIT_THREADS);

        val density:Density!;

        val bas_loc:BasisFunctions!;
        val mol_loc:Molecule[QMAtom]!;

        public def this(mol:Molecule[QMAtom], basisName:String, den:Density) {
            // Console.OUT.println("\tStart making local Molecule and BasisFunctions @ " + here);

            val mol_loc = new Molecule[QMAtom]();
            val nAtoms  = at(mol) { mol.getNumberOfAtoms() };
 
            for((i) in 0..(nAtoms-1)) {
                val sym  = at(mol) { mol.getAtom(i).symbol };
                val centre = at(mol) {mol.getAtom(i).centre };
                mol_loc.addAtom(new QMAtom(sym, centre));
            }

            this.mol_loc = mol_loc;

            val bas_loc = new BasisFunctions(mol_loc, basisName, "basis");
            this.bas_loc = bas_loc;

            // Console.OUT.println("\tDone making local Molecule and BasisFunctions @ " + here);

            // Console.OUT.println("\tMake local copy of density @ " + here);

            val den_loc = new Density(den);
            this.density = den_loc;

            // Console.OUT.println("\tDone making local copy of density @ " + here);

            val maxam = bas_loc.getShellList().getMaximumAngularMomentum();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
                computeInst(i) = new ComputeNew(new TwoElectronIntegrals(maxam), bas_loc.getShellList(), den_loc);          
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
            var jM:Matrix! = new Matrix(N);

            jM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               jM = jM.add(computeInst(i).getJMat());
            } //  end for

            return jM;             
        }

        public def getKMat() : Matrix! {
            val N = density.getRowCount();
            var kM:Matrix! = new Matrix(N);

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

        public def this(bfs:BasisFunctions!, den:Density!) { 
            this.density = den;

            this.shellList = bfs.getShellList();

            val maxam = bfs.getShellList().getMaximumAngularMomentum();
            this.twoEI = new TwoElectronIntegrals(maxam);

            val N = den.getRowCount();

            jM = new Matrix(N);
            kM = new Matrix(N);
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

