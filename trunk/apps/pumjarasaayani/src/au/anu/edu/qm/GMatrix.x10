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
    private global val gMatType : Int;
    private global val computeInst : DistArray[ComputePlace](1){rect,self.at(here)};

    private val bfs : BasisFunctions!;
    private val mol : Molecule[QMAtom]!;

    public def this(N:Int, bfs:BasisFunctions!, mol:Molecule[QMAtom]!, gMatType:Int) {
        super(N);
        this.bfs = bfs;
        this.mol = mol;
        this.gMatType = gMatType;
        val basisName = bfs.getBasisName();

        switch(gMatType) {
        case 0:
        case 1:
            // single-place implementations
            computeInst = null;
            break;
        case 2:
        case 3:
            computeInst = DistArray.make[ComputePlace](Dist.makeUnique(Place.places), ((p) : Point) => new ComputePlaceDirect(N, mol, basisName));
            break;
        case 4:
            computeInst = DistArray.make[ComputePlace](Dist.makeUnique(Place.places), ((p) : Point) => new ComputePlaceFuture(N, mol, basisName));
            break;
        case 5:
        default:
            computeInst = DistArray.make[ComputePlace](Dist.makeUnique(Place.places), ((p) : Point) => new ComputePlaceDirect(N, mol, basisName));
            break;
        } // end switch .. case
    }

    /** top level method to form the G Matrix, depending on gMatType appropriate functions are called */
    public def compute(density:Density!) {
       val timer = new Timer(1);

       timer.start(0);
       switch(gMatType) {
       case 0:
           Console.OUT.println("   GMatrix.computeDirectSerial: " + gMatType);
           computeDirectSerial(density); 
           break;
       case 1:
           Console.OUT.println("   GMatrix.computeDirectLowMemNoAtomic: " + gMatType);
           computeDirectLowMemNoAtomic(density); 
           break;
       case 2:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceNoAtomic: " + gMatType);
           computeDirectMultiPlaceNoAtomic(density);
           break;
       case 3:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceStatic: " + gMatType);
           computeDirectMultiPlaceStatic(density);
           break;
       case 4:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceFuture: " + gMatType);
           computeDirectMultiPlaceFuture(density);
           break;
       case 5:
       default:
           Console.OUT.println("   GMatrix.computeDirectMultiPlaceShellLoop: " + gMatType);
           computeDirectMultiPlaceShellLoop(density);
           break;
       } // end switch .. case
       timer.stop(0);
       Console.OUT.println ("    Time to construct GMatrix: " + (timer.total(0) as Double) / 1e9 + " seconds");
    }

    private def computeDirectSerial(density:Density!) {
        val N = getRowCount();

        makeZero();

        val jMat = new Matrix(N);
        val kMat = new Matrix(N);

        jMat.makeZero();
        kMat.makeZero();

        val shellList = bfs.getShellList();
        val maxam = shellList.getMaximumAngularMomentum();
        val twoE = new TwoElectronIntegrals(maxam);

        val noOfAtoms = mol.getNumberOfAtoms();
        // center a
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
                val iaFunc = aFunc.get(i);

                // center b
                for(var b:Int=0; b<=a; b++) {
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);

                        // center c
                        for(var c:Int=0; c<noOfAtoms; c++) {
                            val cFunc = mol.getAtom(c).getBasisFunctions();
                            val ncFunc = cFunc.size();
                            // basis functions on c
                            for(var k:Int=0; k<ncFunc; k++) {
                                val kcFunc = cFunc.get(k);

                                // center d
                                for(var d:Int=0; d<=c; d++) {
                                    val dFunc = mol.getAtom(d).getBasisFunctions();
                                    val ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(var l:Int=0; l<ndFunc; l++) {
                                        val ldFunc = dFunc.get(l);

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
        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();
        finish foreach((x,y) in gMatrix.region)
                  gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    private def computeDirectLowMemNoAtomic(density:Density!) : void {
        val N = getRowCount();

        makeZero();

        val shellList = bfs.getShellList();

        val computeThreads = Rail.make[ComputeThread!](Runtime.INIT_THREADS);

        val maxam = shellList.getMaximumAngularMomentum();
        for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
            computeThreads(i) = new ComputeThread(N, new TwoElectronIntegrals(maxam), shellList);
        } // end for

        val noOfAtoms = mol.getNumberOfAtoms();
        finish {
          // center a
          for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
               val iaFunc = aFunc.get(i);

               // center b
               for(var b:Int=0; b<=a; b++) {
                   val bFunc = mol.getAtom(b).getBasisFunctions();
                   val nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(var j:Int=0; j<nbFunc; j++) {
                       val jbFunc = bFunc.get(j);

                       // center c
                       for(var c:Int=0; c<noOfAtoms; c++) {
                           val cFunc = mol.getAtom(c).getBasisFunctions();
                           val ncFunc = cFunc.size();
                           // basis functions on c
                           for(var k:Int=0; k<ncFunc; k++) {
                               val kcFunc = cFunc.get(k);

                               // center d
                               for(var d:Int=0; d<=c; d++) {
                                   val dFunc = mol.getAtom(d).getBasisFunctions();
                                   val ndFunc = (d<c) ? dFunc.size() : k+1;
                                   // basis functions on d
                                   for(var l:Int=0; l<ndFunc; l++) {
                                       val ldFunc = dFunc.get(l);

                                       var setIt:Boolean = false;
                                       
                                       outer: while(!setIt) {
                                           for(var idx:Int=0; idx<Runtime.INIT_THREADS; idx++) { 
                                               setIt = computeThreads(idx).setValue(iaFunc, jbFunc, kcFunc, ldFunc);                                               
                                               if (setIt) {
                                                  val idx_loc = idx;
                                                  async computeThreads(idx_loc).compute(density);
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
             val jMatrix = computeThreads(idx).getJMat().getMatrix();
             val kMatrix = computeThreads(idx).getKMat().getMatrix();
             
             finish foreach((x,y) in gMatrix) {
                   gMatrix(x,y) += jMatrix(x,y) - (0.25*kMatrix(x,y));     
             } // finish
        } // end for
    }

    /**
     * Gathers the G matrix contribution from each place
     * and reduces to this place (the complete G matrix).
     */
    public def gatherAndReduceGMatrix() {
        makeZero();
        val N = getRowCount();
        // form the G matrix
        // TODO following need to change once XTENLANG-787 is resolved
        val gMatrix = getMatrix();
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
    }

    private def computeDirectMultiPlaceNoAtomic(density:Density!) {
        val noOfAtoms = mol.getNumberOfAtoms();

        val timer = new Timer(2);

        timer.start(0);

        finish ateach (p in computeInst) {
            computeInst(p).reset(density);
        }

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
                                                    setIt = at(computeInst.dist(p)) { (computeInst(p) as ComputePlaceDirect).setValue(a, b, c, d, i, j, k, l) };

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

        finish ateach (p in computeInst) {
            val comp_loc = computeInst(p) as ComputePlaceDirect!;
            comp_loc.computeGMatContribution();
        }

        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        gatherAndReduceGMatrix();
        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds"); 
    }

    private def computeDirectMultiPlaceStatic(density:Density!) {
        val timer = new Timer(2);

        timer.start(0);

        val noOfAtoms = mol.getNumberOfAtoms();
        val nPlaces = Place.places.length;
        val workPerPlace = noOfAtoms / nPlaces;
        val remainder = noOfAtoms % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);
        Console.OUT.println("\tWork units per place: " + workPerPlace + " remainder " + remainder);

        finish ateach ((placeId) in computeInst) {
            val comp_loc = computeInst(placeId) as ComputePlaceDirect!;
            comp_loc.reset(density);
            val start = placeId < remainder ? (placeId * (workPerPlace + 1)) : (firstChunk + (placeId-remainder) * workPerPlace);
            val end = start + workPerPlace + (placeId < remainder ? 1 : 0);
            comp_loc.compute(start, end);
        } // finish

        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        gatherAndReduceGMatrix();
        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds"); 
    }


    // TODO: need to use shared variable instead
    private global val G = Rail.make[Int](1, (Int)=>0);

    /** Code snippet 3, Bernholdt paper  */
    private def computeDirectMultiPlaceFuture(density:Density!) {
        // init counter
        G(0) = 0;
  
        val timer = new Timer(2);

        timer.start(0);
    
        // center a
        finish ateach ((placeId) in computeInst) {
            val comp_loc = computeInst(placeId) as ComputePlaceFuture!;
            comp_loc.reset(density);

            var myG:Int = 0;
            var L:Int = 0;

            val mol_loc = comp_loc.mol_loc;

            val F1 = future(G) { 
                       var myG:Int; 
                       atomic myG = G(0)++;
                       return myG;
            };

            myG = F1.force();

            val noOfAtoms = at (this) {mol.getNumberOfAtoms()};
            for(var a:Int=0; a<noOfAtoms; a++) {
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
        } // ateach
        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds");
    
        timer.start(1);
        gatherAndReduceGMatrix();
        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds");
    }

    private def computeDirectMultiPlaceShellLoop(density:Density!) {
        val shellList = bfs.getShellList();
        val shellPairs = shellList.getShellPairs();

        val timer = new Timer(2);

        val nPairs = shellPairs.length();
        val nPlaces = Place.places.length;
        val workPerPlace = nPairs / nPlaces;
        val remainder = nPairs % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);
        Console.OUT.println("\tWork units per place: " + workPerPlace + " remainder " + remainder);

        timer.start(0);

        finish ateach ((placeId) in computeInst) {
            val comp_loc = computeInst(placeId) as ComputePlaceDirect!;
            comp_loc.reset(density);
            val start = placeId < remainder ? (placeId * (workPerPlace + 1)) : (firstChunk + (placeId-remainder) * workPerPlace);
            val end = start + workPerPlace + (placeId < remainder ? 1 : 0);
            comp_loc.computeShells(start, end, nPairs);
        }

        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds");

        timer.start(1);
        gatherAndReduceGMatrix();
        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds");

    }

    /** Compute class for the new code - multi place version */
    abstract class ComputePlace {
        var density:Density!;
        val mol_loc:Molecule[QMAtom]!;
        val bas_loc:BasisFunctions!;

        public def this(mol:Molecule[QMAtom], basisName:String) {
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
        }

        /**
         * Prepares to calculate G Matrix using the provided
         * density matrix.  Sets partial contribution to G
         * and other temporary arrays to zero.
         */
        public def reset(density : Density) {
            val den_loc = new Density(density);
            this.density = den_loc;
        }

        public abstract def getGMat() : Matrix!;

        // TODO following method should not be necessary once XTENLANG-787 is resolved
        public def getGMatVal() {
           return getGMat().getValRail();
        }
    }

    class ComputePlaceDirect extends ComputePlace {
        val gMatrixContribution : Matrix!;
        val computeThreads = Rail.make[ComputeThread!](Runtime.INIT_THREADS);

        public def this(N : Int, mol:Molecule[QMAtom], basisName:String) {
            super(mol, basisName);

            gMatrixContribution = new Matrix(N);
        }

        public def reset(density : Density) {
            super.reset(density);
            val maxam = bas_loc.getShellList().getMaximumAngularMomentum();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
                computeThreads(i) = new ComputeThread(gMatrixContribution.getRowCount(), new TwoElectronIntegrals(maxam), bas_loc.getShellList());
            }
        }

        public def getGMat() : Matrix! {
            return gMatrixContribution;
        }

        public def setValue(a:Int, b:Int, c:Int, d:Int, i:Int, j:Int, k:Int, l:Int) : Boolean {
            val iFunc = mol_loc.getAtom(a).getBasisFunctions().get(i);
            val jFunc = mol_loc.getAtom(b).getBasisFunctions().get(j);
            val kFunc = mol_loc.getAtom(c).getBasisFunctions().get(k);
            val lFunc = mol_loc.getAtom(d).getBasisFunctions().get(l);

            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
               val setIt = computeThreads(ix).setValue(iFunc, jFunc, kFunc, lFunc);

               if (setIt) {
                  val ix_loc = ix;
                  async computeThreads(ix_loc).compute(density);
                  return true;
               } // end if
            } // end for

            return false;
        }

        public def compute(i:ContractedGaussian!, j:ContractedGaussian!, 
                           k:ContractedGaussian!, l:ContractedGaussian!) {
            // TODO: this is actually handling only one thread
            for(var ix:Int=0; ix<Runtime.INIT_THREADS; ix++) {
                await (computeThreads(ix).computing == false);

                computeThreads(ix).setValue(i, j, k, l);
                val ix_loc = ix;
                async computeThreads(ix_loc).compute(density);
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
                                       computeThreads(0).computeSingle(iaFunc, jbFunc, kcFunc, ldFunc, density);
                                   } // end l
                               } // center d
                           } // end k
                       } // center c
                   } // end j
               } // center b
             } // end i
          } // center a

            computeGMatContribution();
        }

        public def computeShells(startShell:Int, endShell:Int, nPairs:Int) {
            val bfs = computeThreads(0).shellList.getShellPrimitives();
            val shellPairs = computeThreads(0).shellList.getShellPairs();

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
                    computeThreads(0).computeSingle2(aFunc, bFunc, cFunc, dFunc,
                                radiusABSquared, aAng, bAng, cAng, dAng, angMomAB,
                                aStrt, bStrt, cStrt, dStrt, aLim, bLim, abLim,
                                density);
                } // end for
            } // end for

            computeGMatContribution();
        }


        private def getJMat() : Matrix! {
            val N = density.getRowCount();
            var jM:Matrix! = new Matrix(N);

            jM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               jM = jM.add(computeThreads(i).getJMat());
            } //  end for

            return jM;             
        }

        private def getKMat() : Matrix! {
            val N = density.getRowCount();
            var kM:Matrix! = new Matrix(N);

            kM.makeZero();
            for(var i:Int=0; i<Runtime.INIT_THREADS; i++) {
               kM = kM.add(computeThreads(i).getKMat());
            } //  end for

            return kM;
        }

        /**
         * Computes this place's contribution to the G Matrix
         * given previously calculated J, K matrices.
         */
        public def computeGMatContribution() {
            val jMat = getJMat().getMatrix();
            val kMat = getKMat().getMatrix();
            val gMat = gMatrixContribution.getMatrix();
            for (p in jMat.region) {
                gMat(p) = jMat(p) - 0.25 * kMat(p);
            }
        }
    }

    class ComputePlaceFuture extends ComputePlace {
        val twoEI:TwoElectronIntegrals!;
        val jM:Matrix!;
        val kM:Matrix!;

        public def this(N : Int, mol:Molecule[QMAtom], basisName:String) {
            super(mol, basisName);

            val maxam = bas_loc.getShellList().getMaximumAngularMomentum();
            this.twoEI = new TwoElectronIntegrals(maxam);

            jM = new Matrix(N);
            kM = new Matrix(N); 
         }

        public def reset(density : Density) {
            super.reset(density);
            jM.makeZero();
            kM.makeZero();
        }

         public def compute2EAndRecord(iaFunc:ContractedGaussian!, jbFunc:ContractedGaussian!, 
                                       kcFunc:ContractedGaussian!, ldFunc:ContractedGaussian!) {
             twoEI.compute2EAndRecord(iaFunc, jbFunc, kcFunc, ldFunc, bas_loc.getShellList(), jM, kM, density);
         }

        public def getGMat() : Matrix! {
            val gM = new Matrix(jM.getRowCount());
            val gMat = gM.getMatrix();
            val jMat = jM.getMatrix();
            val kMat = kM.getMatrix();

            for (p in jMat.region) {
                gMat(p) = jMat(p) - 0.25 * kMat(p);
            }
            return gM;
        }
    }

    class ComputeThread {
        var computing:Boolean = false;

        // row count for density and G matrices
        val N : Int;

        var i:ContractedGaussian!, j:ContractedGaussian!, 
            k:ContractedGaussian!, l:ContractedGaussian!;

        val twoEI:TwoElectronIntegrals!;
        val shellList:ShellList{self.at(this)};
        val jMat:Matrix!, kMat:Matrix!;

        public def this(N: Int, te:TwoElectronIntegrals!, sh:ShellList) {
            this.N = N;
            this.twoEI = te;
            this.shellList = sh;

            jMat = new Matrix(N);
            kMat = new Matrix(N);

            jMat.makeZero();
            kMat.makeZero();
        }

        public def compute(density : Density!) {
            twoEI.compute2EAndRecord(i, j, 
                                     k, l,  
                                     shellList, 
                                     jMat, kMat, density);
            atomic computing = false;
        }

        public def computeSingle(i:ContractedGaussian!, j:ContractedGaussian!,
                                 k:ContractedGaussian!, l:ContractedGaussian!,
                                 density : Density!) {
            twoEI.compute2EAndRecord(i, j, k, l, 
                                     shellList,
                                     jMat, kMat, density);
        }

        public def computeSingle2(i:ContractedGaussian!, j:ContractedGaussian!,
                                  k:ContractedGaussian!, l:ContractedGaussian!, 
                                  radiusABSquared:Double,
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int, abLim:Int,
                                  density : Density!) {
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
}

