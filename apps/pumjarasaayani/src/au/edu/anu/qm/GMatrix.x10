/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.util.ArrayList;
import x10.util.Team;
import x10.util.concurrent.Future;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;

/**
 * GMatrix.x10
 *
 * GMatrix in HF calculation
 *
 * @author: V.Ganesh
 */
public class GMatrix extends Matrix {
    private val gMatType : Int;
    private val computeInst : DistArray[ComputePlace](1);

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom], gMatType:Int) {
        super(N);
        this.bfs = bfs;
        this.mol = molecule;
        this.gMatType = gMatType;

        val basisName = bfs.getBasisName();

        computeInst = DistArray.make[ComputePlace](Dist.makeUnique(), (Point) => new ComputePlace(N, molecule, bfs));
    }

    /** top level method to form the G Matrix, depending on gMatType appropriate functions are called */
    public def compute(density:Density) {
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

    private def computeDirectSerial(density:Density) {
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
        // centre a
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
                val iaFunc = aFunc.get(i);

                // centre b
                for(var b:Int=0; b<=a; b++) {
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);

                        // centre c
                        for(var c:Int=0; c<noOfAtoms; c++) {
                            val cFunc = mol.getAtom(c).getBasisFunctions();
                            val ncFunc = cFunc.size();
                            // basis functions on c
                            for(var k:Int=0; k<ncFunc; k++) {
                                val kcFunc = cFunc.get(k);

                                // centre d
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
                                } // centre d
                            } // end k
                        } // centre c
                    } // end j
                } // centre b
            } // end i
        } // centre a
     
        // form the G matrix
        val gMatrix = getMatrix();
        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();
        for([x,y] in gMatrix.region) gMatrix(x,y) = jMatrix(x,y) - (0.25*kMatrix(x,y));     
    }

    private def computeDirectLowMemNoAtomic(density:Density) : void {
        val N = getRowCount();

        makeZero();

        val shellList = bfs.getShellList();

        val computeThreads = new Array[ComputeThread](Runtime.NTHREADS);

        for(var i:Int=0; i<Runtime.NTHREADS; i++) {
            computeThreads(i) = new ComputeThread(N, shellList);
        } // end for

        val noOfAtoms = mol.getNumberOfAtoms();
        finish {
          // centre a
          for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
               val iaFunc = aFunc.get(i);

               // centre b
               for(var b:Int=0; b<=a; b++) {
                   val bFunc = mol.getAtom(b).getBasisFunctions();
                   val nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(var j:Int=0; j<nbFunc; j++) {
                       val jbFunc = bFunc.get(j);

                       // centre c
                       for(var c:Int=0; c<noOfAtoms; c++) {
                           val cFunc = mol.getAtom(c).getBasisFunctions();
                           val ncFunc = cFunc.size();
                           // basis functions on c
                           for(var k:Int=0; k<ncFunc; k++) {
                               val kcFunc = cFunc.get(k);

                               // centre d
                               for(var d:Int=0; d<=c; d++) {
                                   val dFunc = mol.getAtom(d).getBasisFunctions();
                                   val ndFunc = (d<c) ? dFunc.size() : k+1;
                                   // basis functions on d
                                   for(var l:Int=0; l<ndFunc; l++) {
                                       val ldFunc = dFunc.get(l);

                                       var setIt:Boolean = false;
                                       
                                       outer: while(!setIt) {
                                           for(var idx:Int=0; idx<Runtime.NTHREADS; idx++) { 
                                               setIt = computeThreads(idx).setValue(iaFunc, jbFunc, kcFunc, ldFunc);                                               
                                               if (setIt) {
                                                  val idx_loc = idx;
                                                  async computeThreads(idx_loc).compute(density);
                                                  break outer;
                                               } // end if
                                           } // end for
                                       } // end while
				   } // end l
			       } // centre d
                           } // end k
                       } // centre c
                   } // end j
               } // centre b
            } // end i
          } // centre a
        } // finish
        
        // form the G matrix
        val gMatrix = getMatrix();
        for(var idx:Int=0; idx<Runtime.NTHREADS; idx++) { 
             val jMatrix = computeThreads(idx).getJMat().getMatrix();
             val kMatrix = computeThreads(idx).getKMat().getMatrix();
             
             for([x,y] in gMatrix) {
                   gMatrix(x,y) += jMatrix(x,y) - (0.25*kMatrix(x,y));     
             }
        } // end for
    }

    private def computeDirectMultiPlaceNoAtomic(density:Density) {
        val noOfAtoms = mol.getNumberOfAtoms();

        val timer = new Timer(2);

        timer.start(0);

        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish ateach (p in computeInst) {
            computeInst(p).reset(density);
        }

        finish {
            // centre a
            for (a in 0..(noOfAtoms-1)) {
                val aFunc = mol.getAtom(a).getBasisFunctions();
                val naFunc = aFunc.size();
                // basis functions on a
                for (i in 0..(naFunc-1)) {
                    val iaFunc = aFunc.get(i);

                    // centre b
                    for (b in 0..a) {
                        val bFunc = mol.getAtom(b).getBasisFunctions();
                        val nbFunc = (b<a) ? bFunc.size() : i+1;
                        // basis functions on b
                        for (j in 0..(nbFunc-1)) {
                            val jbFunc = bFunc.get(j);

                            // centre c
                            for (c in 0..(noOfAtoms-1)) {
                                val cFunc = mol.getAtom(c).getBasisFunctions();
                                val ncFunc = cFunc.size();
                                // basis functions on c
                                for (k in 0..(ncFunc-1)) {
                                    val kcFunc = cFunc.get(k);

                                    // centre d
                                    for (d in 0..c) {
                                       val dFunc = mol.getAtom(d).getBasisFunctions();
                                       val ndFunc = (d<c) ? dFunc.size() : k+1;
                                        // basis functions on d
                                        for (l in 0..(ndFunc-1)) {
                                            val ldFunc = dFunc.get(l);
                                            var setIt:Boolean = false;

                                            outer: while(!setIt) {
                                                for(p in computeInst) {
                                                    setIt = at(computeInst.dist(p)) { computeInst(p).setValue(a, b, c, d, i, j, k, l) };

                                                    if (setIt) break outer;
                                                } // end for
                                            } // end while                               
                                        } // end l
                                    } // centre d
                                } // end k
                            } // centre c
                        } // end j
                    } // centre b
                } // end i
            } // centre a
        } // finish

        timer.stop(0);
        Console.OUT.println("\tTime for actual computation: " + (timer.total(0) as Double) / 1e9 + " seconds"); 

        timer.start(1);
        gatherAndReduceGMatrix();
        timer.stop(1);
        Console.OUT.println("\tTime for summing up GMatrix bits: " + (timer.total(1) as Double) / 1e9 + " seconds"); 
    }

    /**
     * Gathers the G matrix contribution from each place
     * and reduces to this place (the complete G matrix).
     */
    public def gatherAndReduceGMatrix() {
        makeZero();
        val N = getRowCount();

        val sum = (a:Double, b:Double) => (a+b);
        // form the G matrix
        val gMatrix = getMatrix();
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for(p in computeInst) async {
            val gVal = at(computeInst.dist(p)) { computeInst(p).getGMatContributionArray() };
            atomic { gMatrix.map[Double,Double](gMatrix, gVal, sum); }
        } // end for
    }

    private def computeDirectMultiPlaceStatic(density:Density) {
        val noOfAtoms = mol.getNumberOfAtoms();
        val nPlaces = Place.MAX_PLACES;
        val workPerPlace = noOfAtoms / nPlaces;
        val remainder = noOfAtoms % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);
        Console.OUT.println("\tWork units per place: " + workPerPlace + " remainder " + remainder);

        val gMat = getMatrix();
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at (Place.place(placeId)) {
                //val placeTimer = new Timer(1);
                //placeTimer.start(0);

                val comp_loc = computeInst(placeId);
                comp_loc.reset(density);
                val start = placeId < remainder ? (placeId * (workPerPlace + 1)) : (firstChunk + (placeId-remainder) * workPerPlace);
                val end = start + workPerPlace + (placeId < remainder ? 1 : 0);
                comp_loc.compute(start, end);

                //placeTimer.stop(0);
                //Console.OUT.println("\tcompute at " + here + " " + (placeTimer.total(0) as Double) / 1e9 + " seconds");

                comp_loc.getGMatContributionArray()
            };

            // gather and reduce my gMatrix contribution
            //val gatherTimer = new Timer(1);
            //gatherTimer.start(0);
            val sum = (a:Double, b:Double) => (a+b);
            atomic { gMat.map[Double,Double](gMat, placeContribution, sum); }
            //gatherTimer.stop(0);
            //Console.OUT.println("\tgather from " + placeId + " " + (gatherTimer.total(0) as Double) / 1e9 + " seconds");
        }
    }

    /** Code snippet 3, Bernholdt paper  */
    private def computeDirectMultiPlaceFuture(density:Density) {
        val G = new SharedCounter();
        G.set(Place.MAX_PLACES); // each place is first assigned the counter value of its own place number
        makeZero();
        val gMat = getMatrix();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at (Place.place(placeId)) {
                val placeTimer = new Timer(1);
                placeTimer.start(0);

                val comp_loc = computeInst(placeId);
                comp_loc.reset();

                var myG : Int = placeId;

                val bfs = comp_loc.computeThreads(0).shellList.getShellPrimitives();
                val nPrimitives = comp_loc.computeThreads(0).shellList.getNumberOfShellPrimitives();
                val nPairs = comp_loc.computeThreads(0).shellList.getNumberOfShellPairs();

                var completedHere : Int = 0;
                for(var i:Int=myG; i<nPairs; i++) {
                    if (i == myG) {
                        val F2 = Future.make[Int](() => G.getAndIncrement());
                        comp_loc.computeOneShellPair(i, nPrimitives, bfs);
                        completedHere++;
                        myG = F2.force();
                    }
                }

                placeTimer.stop(0);
                Console.OUT.println("\tcompute at " + here + " completed " + completedHere + " " + (placeTimer.total(0) as Double) / 1e9 + " seconds");

                comp_loc.getGMatContributionArray()
            };

            // gather and reduce my gMatrix contribution
            //val gatherTimer = new Timer(1);
            //gatherTimer.start(0);
            val sum = (a:Double, b:Double) => (a+b);
            atomic { gMat.map[Double,Double](gMat, placeContribution, sum); }
            //gatherTimer.stop(0);
            //Console.OUT.println("\tgather from " + placeId + " " + (gatherTimer.total(0) as Double) / 1e9 + " seconds");
        } // ateach
    }

    private def computeDirectMultiPlaceShellLoop(density:Density) {
        val shellList = bfs.getShellList();
        val nPairs = shellList.getNumberOfShellPairs();
        val nPlaces = Place.MAX_PLACES;
        val workPerPlace = nPairs / nPlaces;
        val remainder = nPairs % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);
        Console.OUT.println("\tWork units per place: " + workPerPlace + " remainder " + remainder);

        makeZero();
        val gMat = getMatrix();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at (Place.place(placeId)) {
                //val placeTimer = new Timer(1);
                //placeTimer.start(0);
                val comp_loc = computeInst(placeId);
                comp_loc.reset();
                comp_loc.computeShells(nPairs);
                //placeTimer.stop(0);
                //Console.OUT.println("\tcompute at " + here + " " + (placeTimer.total(0) as Double) / 1e9 + " seconds");
                comp_loc.getGMatContributionArray()
            };

            // gather and reduce my gMatrix contribution
            //val gatherTimer = new Timer(1);
            //gatherTimer.start(0);
            val sum = (a:Double, b:Double) => (a+b);
            atomic { gMat.map[Double,Double](gMat, placeContribution, sum); }
            //gatherTimer.stop(0);
            //Console.OUT.println("\tgather at " + here + " " + (gatherTimer.total(0) as Double) / 1e9 + " seconds");
        }
    }

    /** Compute class for the new code - multi place version */
    static class ComputePlace {
        var density:Density;
        val mol_loc:Molecule[QMAtom];

        val gMatrixContribution : Matrix;
        val jMatrixContribution : Matrix;
        val kMatrixContribution : Matrix;

        public val computeThreads = new Array[ComputeThread](Runtime.NTHREADS);

        public def this(N : Int, mol:Molecule[QMAtom], bfs:BasisFunctions) {
            this.mol_loc = mol;

            gMatrixContribution = new Matrix(N);
            jMatrixContribution = new Matrix(N);
            kMatrixContribution = new Matrix(N);

            // TODO not required with broadcast ateach above - XTENLANG-1725
            val noOfElectrons = mol.getNumberOfElectrons();
            val noOfOccupancies = noOfElectrons / 2;
            density = new Density(N, noOfOccupancies);

            val shellList = bfs.getShellList();
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i) = new ComputeThread(gMatrixContribution.getRowCount(), shellList);
            }
        }

        /**
         * Prepares to calculate G Matrix using the provided
         * density matrix.  Sets partial contribution to G
         * and other temporary arrays to zero.
         */
        public def reset(density : Density) {
            this.density = density;
            gMatrixContribution.makeZero();
            jMatrixContribution.makeZero();
            kMatrixContribution.makeZero();
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i).reset();
            }
        }

        /**
         * Prepares to calculate G Matrix using a broadcast
         * of the density matrix from place 0.  Sets partial 
         * contribution to G and other temporary arrays to zero.
         */
        public def reset() {
            // TODO should use ateach broadcast - XTENLANG-1725
            val densityMatrix = density.getMatrix();
            Team.WORLD.bcast[Double](here.id, 0, densityMatrix, 0, densityMatrix, 0, densityMatrix.size);
            gMatrixContribution.makeZero();
            jMatrixContribution.makeZero();
            kMatrixContribution.makeZero();
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i).reset();
            }
        }

        /**
         * Computes this place's contribution to the G Matrix
         * given previously calculated J, K matrices.
         */
        public def getGMatContributionArray() {
            computeJMatContribution();
            computeKMatContribution();

            val jMat = jMatrixContribution.getMatrix();
            val kMat = kMatrixContribution.getMatrix();
            val gMat = gMatrixContribution.getMatrix();
            for (p in jMat.region) {
                gMat(p) = jMat(p) - 0.25 * kMat(p);
            }
           return gMatrixContribution.getMatrix();
        }

        public def computeOneShellPair(shellPrim : Int, nPrimitives : Int, bfs : ArrayList[ContractedGaussian]) {
            val a = shellPrim / nPrimitives;
            val b = shellPrim % nPrimitives;

            val aFunc = bfs(a);
            val bFunc = bfs(b);

            val aStrt = aFunc.intIndex;
            val bStrt = bFunc.intIndex;
            val aAng  = aFunc.getMaximumAngularMomentum();
            val bAng  = bFunc.getMaximumAngularMomentum();

            val aa = aStrt + aAng;
            val bb = bStrt + bAng;

            if (aa < bb) return;

            val angMomAB = aAng + bAng;
            val aLim = ((aAng+1)*(aAng+2)/2);
            val bLim = ((bAng+1)*(bAng+2)/2);
            val abLim = aLim * bLim;

            val radiusABSquared = aFunc.distanceSquaredFrom(bFunc);

            for (c in 0..(nPrimitives-1)) {
                val cFunc = bfs(c);
                val cStrt = cFunc.intIndex;
                val cAng  = cFunc.getMaximumAngularMomentum();
                val cc = cStrt + cAng;

                if (aa >= cc) {
                    for (d in 0..(nPrimitives-1)) {
                        val dFunc = bfs(d);
                        val dStrt = dFunc.intIndex;
                        val dAng  = dFunc.getMaximumAngularMomentum();
                        val dd = dStrt + dAng;

                        if (cc >= dd) {
                            computeThreads(0).computeSingle2(aFunc, bFunc, cFunc, dFunc,
                                        radiusABSquared, aAng, bAng, cAng, dAng, angMomAB,
                                        aStrt, bStrt, cStrt, dStrt, aLim, bLim, abLim,
                                        density);
                        }
                    }
                    if (here == Place.FIRST_PLACE) {
                        // hack - first place holds shared counter
                        // check for incoming messages (e.g. shared counter updates)
                        // TODO instead, split work into smaller activities
                        Runtime.probe();
                    }
                }
            }
        }

        public def setValue(a:Int, b:Int, c:Int, d:Int, i:Int, j:Int, k:Int, l:Int) : Boolean {
            val iFunc = mol_loc.getAtom(a).getBasisFunctions().get(i);
            val jFunc = mol_loc.getAtom(b).getBasisFunctions().get(j);
            val kFunc = mol_loc.getAtom(c).getBasisFunctions().get(k);
            val lFunc = mol_loc.getAtom(d).getBasisFunctions().get(l);

            val computeThreads = this.computeThreads; // TODO this should not be required XTENLANG-1913
            for(var ix:Int=0; ix<Runtime.NTHREADS; ix++) {
               val setIt = computeThreads(ix).setValue(iFunc, jFunc, kFunc, lFunc);

               if (setIt) {
                  val ix_loc = ix;
                  async computeThreads(ix_loc).compute(density);
                  return true;
               } // end if
            } // end for

            return false;
        }

        public def compute(start:Int, end:Int) {
            val noOfAtoms = mol_loc.getNumberOfAtoms();

            for(var a:Int=start; a<end; a++) {
             val aFunc = mol_loc.getAtom(a).getBasisFunctions();
             val naFunc = aFunc.size();
             // basis functions on a
             for(var i:Int=0; i<naFunc; i++) {
               val iaFunc = aFunc.get(i);

               // centre b
               for(var b:Int=0; b<=a; b++) {
                   val bFunc = mol_loc.getAtom(b).getBasisFunctions();
                   val nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(var j:Int=0; j<nbFunc; j++) {
                       val jbFunc = bFunc.get(j);

                       // centre c
                       for(var c:Int=0; c<noOfAtoms; c++) {
                           val cFunc = mol_loc.getAtom(c).getBasisFunctions();
                           val ncFunc = cFunc.size();
                           // basis functions on c
                           for(var k:Int=0; k<ncFunc; k++) {
                               val kcFunc = cFunc.get(k);

                               // centre d
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
                               } // centre d
                           } // end k
                       } // centre c
                   } // end j
               } // centre b
             } // end i
          } // centre a
        }

        public def computeShells(nPairs:Int) {
            val bfs = computeThreads(0).shellList.getShellPrimitives();
            val nPrimitives = computeThreads(0).shellList.getNumberOfShellPrimitives();

            for(var i:Int=here.id; i<nPairs; i+=Place.MAX_PLACES) {
                computeOneShellPair(i, nPrimitives, bfs);
            }
        }

        private def computeJMatContribution() : Matrix {
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
               jMatrixContribution.addInPlace(computeThreads(i).getJMat());
            }

            return jMatrixContribution;             
        }

        private def computeKMatContribution() : Matrix {
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
               kMatrixContribution.addInPlace(computeThreads(i).getKMat());
            }

            return kMatrixContribution;
        }
    }

    static class ComputeThread {
        var computing:Boolean = false;

        // row count for density and G matrices
        val N : Int;

        var i:ContractedGaussian, j:ContractedGaussian, 
            k:ContractedGaussian, l:ContractedGaussian;

        val twoEI:TwoElectronIntegrals;
        val shellList:ShellList;
        val jMat:Matrix, kMat:Matrix;

        public def this(N: Int, sh:ShellList) {
            this.N = N;
            this.twoEI = new TwoElectronIntegrals(sh.getMaximumAngularMomentum());
            this.shellList = sh;

            jMat = new Matrix(N);
            kMat = new Matrix(N);
        }

        public def reset() {
            jMat.makeZero();
            kMat.makeZero();
        }

        public def compute(density : Density) {
            twoEI.compute2EAndRecord(i, j, 
                                     k, l,  
                                     shellList, 
                                     jMat, kMat, density);
            atomic computing = false;
        }

        public def computeSingle(i:ContractedGaussian, j:ContractedGaussian,
                                 k:ContractedGaussian, l:ContractedGaussian,
                                 density : Density) {
            twoEI.compute2EAndRecord(i, j, k, l, 
                                     shellList,
                                     jMat, kMat, density);
        }

        public def computeSingle2(i:ContractedGaussian, j:ContractedGaussian,
                                  k:ContractedGaussian, l:ContractedGaussian, 
                                  radiusABSquared:Double,
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int, abLim:Int,
                                  density : Density) {
            twoEI.compute2EAndRecord2(i, j, k, l,
                                      shellList,
                                      jMat, kMat, density,
                                      radiusABSquared,
                                      aAng, bAng, cAng, dAng, angMomAB,
                                      aStrt, bStrt, cStrt, dStrt,
                                      aLim, bLim, abLim);
        }
        
        public def setValue(i:ContractedGaussian, j:ContractedGaussian,
                            k:ContractedGaussian, l:ContractedGaussian) : Boolean {
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

