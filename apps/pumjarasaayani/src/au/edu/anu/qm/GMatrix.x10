/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.util.ArrayList;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10.matrix.DenseMatrix;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;

/**
 * G matrix in HF calculation
 * Integral screening by cutoff based on Häser & Ahlrichs (1989)
 * @see Häser, M. and Ahlrichs, R. (1989). "Improvements on the Direct SCF
 *   method". J. Comp. Chem. 10 (1) pp.104-111.
 *
 * @author: V.Ganesh
 */
public class GMatrix extends DenseMatrix{self.M==self.N} {
    public static DEFAULT_GMATTYPE=4;

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val gMatType : Int;
    private val computeInst : DistArray[ComputePlace](1);

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    static THRESH:Double = 1.0e-8;
    private var thresh2:Double = 0.0;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom]):GMatrix{self.M==N,self.N==N} {
        super(N, N);
        this.bfs = bfs;
        this.mol = molecule;

        val jd = JobDefaults.getInstance();
        this.gMatType = jd.gMatrixParallelScheme;

        val nPlaces = Place.MAX_PLACES;
        switch(gMatType) {
        case 0:
            Console.OUT.println("GMatrix.computeSerial: " + gMatType);
            break;
        case 1:
            Console.OUT.println("GMatrix.computeThreadedLowMemNoAtomicByBF: " + gMatType);
            break;
        case 2:
            Console.OUT.println("GMatrix.computeDistNoAtomicByBF: " + gMatType);
            break;
        case 3:
            Console.OUT.println("GMatrix.computeDistStaticByAtoms: " + gMatType);
            val noOfAtoms = mol.getNumberOfAtoms();
            Console.OUT.println("\tAtoms per place: " + (noOfAtoms / nPlaces) + " remainder " + (noOfAtoms % nPlaces));
            break;
        case 4:
            Console.OUT.println("GMatrix.computeDistDynamicByShells: " + gMatType);
            break;
        case 6:
            Console.OUT.println("with allreduce:");
        case 5:
        default:
            Console.OUT.println("GMatrix.computeDistStaticByShells: " + gMatType);
            val shellList = bfs.getShellList();
            val nPairs = shellList.getNumberOfShellPairs();
            Console.OUT.println("\tShell pairs per place: " + (nPairs / nPlaces) + " remainder " + (nPairs % nPlaces));
            break;
        case 7:
            Console.OUT.println("GMatrix.computeDistThreadedDynamicByShells: " + gMatType);
            break;
        }

        val dCut = new DenseMatrix(N, N);
        val qCut = new DenseMatrix(N, N);
        // Schwarz cut-off: Häser & Ahlrichs eqn 12
        val shellList = bfs.getShellList();
        val maxam = shellList.getMaximumAngularMomentum();
        val twoE = new TwoElectronIntegrals(maxam, bfs.getNormalizationFactors(), THRESH);


        val fakeDensity = new Density(N, 2, 1.0);

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
                        twoE.compute2EAndRecord(iaFunc, jbFunc, iaFunc, jbFunc, shellList, qCut, dCut, fakeDensity);
                    }
                }
            }
        }

        var maxEst:Double = 0.0;
        for(var a:Int=0; a<N; a++) for(var b:Int=0; b<N; b++) {
            if (a==b) qCut(a,b)*=.5; else qCut(a,b)*=.25;
            qCut(a,b)=Math.sqrt(Math.abs(qCut(a,b)));
            if (qCut(a,b)>maxEst) maxEst=qCut(a,b);
            //Console.OUT.printf("%d %d %e \n",a,b,qCut(a,b));
        }

        Console.OUT.printf("\tmaxEst %.4g\n", maxEst);
        val maxEstVal = maxEst;

        computeInst = DistArray.make[ComputePlace](Dist.makeUnique(), (Point) => new ComputePlace(N, molecule, bfs, qCut, dCut, maxEstVal));
    }

    /** top level method to form the G Matrix, depending on gMatType appropriate functions are called */
    public def compute(density:Density) {
       timer.start(0);
       switch(gMatType) {
       case 0:
           computeSerial(density); 
           break;
       case 1:
           computeThreadedLowMemNoAtomicByBF(density); 
           break;
       case 2:
           computeDistNoAtomicByBF(density);
           break;
       case 3:
           computeDistStaticByAtoms(density);
           break;
       case 4:
       default:
           computeDistDynamicByShells(density);
           break;
       case 5:
           computeDistStaticByShells(density);
           break;
       case 6:
           computeDistStaticByShellsReduce(density);
           break;
       case 7:
           computeDistThreadedDynamicByShells(density);
           break;
       } // end switch .. case
       timer.stop(0);
       Console.OUT.printf("    Time to construct GMatrix: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }

    private def computeSerial(density:Density) {
        val computePlace = computeInst(here.id);
        computePlace.reset(density);
        val computeThread = computePlace.computeThreads(0);

        val shellList = bfs.getShellList();
        val shellPrimitives = shellList.getShellPrimitives();
        val numPrimitives = shellList.getNumberOfShellPrimitives();

        var totInt:Long=0;

        for(var a:Int=0; a<numPrimitives; a++) {
            val aFunc = shellPrimitives(a);
            for(var b:Int=0; b<numPrimitives; b++) {
                val bFunc = shellPrimitives(b);
                totInt += computeThread.computeOneShellPair(aFunc, bFunc, numPrimitives, shellPrimitives);
            }
        }
        Console.OUT.println("    totInt = "+totInt);
        // form the G matrix
        val jMat = computeThread.getJMat();
        val kMat = computeThread.getKMat();

        for([x,y] in 0..(M-1)*0..(N-1)) {
           this(x,y) = jMat(x,y) - (0.25*kMat(x,y));
        }

        //Console.OUT.println("J Mat");
        //Console.OUT.println(computeThread.getJMat());

        //Console.OUT.println("K Mat");
        //Console.OUT.println(computeThread.getKMat());
// vvvv For development purpose vvvvvv
        val jMatrix=computeThread.getJMat();
        val kMatrix=computeThread.getKMat();

        val eJ = density.clone().mult(density, jMatrix).trace();
        Console.OUT.printf("  EJ = %.6f a.u.\n", eJ);

        val eK = density.clone().mult(density, kMatrix).trace();
        Console.OUT.printf("  EK = %.6f a.u.\n", eK);
// ^^^^ It is not required for normal calculation ^^^^^
    }

    private def computeThreadedLowMemNoAtomicByBF(density:Density) : void {
        reset();

        val computePlace = computeInst(here.id);
        computePlace.reset(density);

        val computeThreads = computePlace.computeThreads;

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
        for(var idx:Int=0; idx<Runtime.NTHREADS; idx++) { 
            val jMat = computeThreads(idx).getJMat();
            val kMat = computeThreads(idx).getKMat();
             
            for([x,y] in 0..(M-1)*0..(N-1)) {
               this(x,y) = jMat(x,y) - (0.25*kMat(x,y));
            }
        } // end for
    }

    private def computeDistNoAtomicByBF(density:Density) {
        val noOfAtoms = mol.getNumberOfAtoms();

        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish ateach(p in computeInst) {
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

        gatherAndReduceGMatrix();
    }

    /**
     * Gathers the G matrix contribution from each place
     * and reduces to this place (the complete G matrix).
     */
    public def gatherAndReduceGMatrix() {
        reset();
        // form the G matrix
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for(p in computeInst) async {
            val gVal = at(computeInst.dist(p)) { computeInst(p).getGMatContributionArray() };
            atomic super.cellAdd(gVal);
        } // end for
    }

    private def computeDistStaticByAtoms(density:Density) {
        reset();
        val noOfAtoms = mol.getNumberOfAtoms();
        val nPlaces = Place.MAX_PLACES;
        val workPerPlace = noOfAtoms / nPlaces;
        val remainder = noOfAtoms % nPlaces;
        val firstChunk = remainder * (workPerPlace + 1);

        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at(Place.place(placeId)) {
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
            atomic super.cellAdd(placeContribution);
            //gatherTimer.stop(0);
            //Console.OUT.println("\tgather from " + placeId + " " + (gatherTimer.total(0) as Double) / 1e9 + " seconds");
        }
    }

    /** Code snippet 3, Bernholdt paper  */
    private def computeDistDynamicByShells(density:Density) {
        val G = new SharedCounter();
        G.set(Place.MAX_PLACES); // each place is first assigned the counter value of its own place number

        reset();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at(Place.place(placeId)) {
                val placeTimer = new Timer(1);
                placeTimer.start(0);

                val comp_loc = computeInst(placeId);
                comp_loc.reset();
                val computeThread = comp_loc.computeThreads(0);

                var myG : Int = placeId;

                val shellList = computeThread.shellList;
                val shellPrimitives = shellList.getShellPrimitives();
                val numPrimitives = shellList.getNumberOfShellPrimitives();
                val numPairs = shellList.getNumberOfShellPairs();

                var placePairs:Int = 0;
                var placeInts:Long = 0;
                for(var i:Int=myG; i<numPairs; i++) {
                    if (i == myG) {
                        finish {
                          async placeInts += computeThread.computeOneShellPair(i, numPrimitives, shellPrimitives);
                          myG = G.getAndIncrement();
                        }
                        placePairs++;
                    }
                }
                placeTimer.stop(0);
                Console.OUT.printf("\tcompute at %s pairs %i integrals %i %.4g seconds\n", here, placePairs, placeInts, ((placeTimer.total(0) as Double) / 1e9));

                comp_loc.getGMatContributionArray()
            };
            // gather and reduce my gMatrix contribution
            //val gatherTimer = new Timer(1);
            //gatherTimer.start(0);
            atomic super.cellAdd(placeContribution);
            //gatherTimer.stop(0);
            //Console.OUT.printf("\tgather from %i %.3g seconds\n", placeId, ((gatherTimer.total(0) as Double) / 1e9));
        }
    }

    private def computeDistStaticByShells(density:Density) {
        val shellList = bfs.getShellList();
        val nPairs = shellList.getNumberOfShellPairs();

        reset();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at(Place.place(placeId)) {
                //val placeTimer = new Timer(1);
                //placeTimer.start(0);
                val comp_loc = computeInst(placeId);
                comp_loc.reset();
                val totInt = comp_loc.computeShells(nPairs);
                //placeTimer.stop(0);
                //Console.OUT.println("    totInt at " + here + " = " + totInt);
                //Console.OUT.println("\tcompute at " + here + " " + (placeTimer.total(0) as Double) / 1e9 + " seconds");
                comp_loc.getGMatContributionArray()
            };

            // gather and reduce my gMatrix contribution
            val gatherTimer = new Timer(1);
            gatherTimer.start(0);
            atomic super.cellAdd(placeContribution);
            gatherTimer.stop(0);
            Console.OUT.printf("\tgather from %i %.3g seconds\n", placeId, ((gatherTimer.total(0) as Double) / 1e9));
        }
    }

    private def computeDistStaticByShellsReduce(density:Density) {
        val shellList = bfs.getShellList();
        val nPairs = shellList.getNumberOfShellPairs();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish ateach ([placeId] in computeInst) {
            //val placeTimer = new Timer(1);
            //placeTimer.start(0);
            val comp_loc = computeInst(placeId);
            comp_loc.reset();
            val totInt = comp_loc.computeShells(nPairs);
            //placeTimer.stop(0);
            //Console.OUT.println("    totInt at " + here + " = " + totInt);
            //Console.OUT.println("\tcompute at " + here + " " + (placeTimer.total(0) as Double) / 1e9 + " seconds");
            comp_loc.allreduceGMat();
        }
        // TODO should reduce directly to gMat
        DenseMatrix.copy(computeInst(0).gMatrixContribution, this);
    }

    /** As per computeDistDynamicByShells with multithreaded places  */
    private def computeDistThreadedDynamicByShells(density:Density) {
        val G = new SharedCounter();
        G.set(Place.MAX_PLACES); // each place is first assigned the counter value of its own place number

        reset();

        computeInst(0).density = density; // prepare for broadcast
        val computeInst = this.computeInst; // TODO this should not be required XTENLANG-1913
        finish for ([placeId] in computeInst) async {
            val placeContribution = at(Place.place(placeId)) {
                val placeTimer = new Timer(1);
                placeTimer.start(0);

                val comp_loc = computeInst(placeId);
                comp_loc.reset();

                val shellList = comp_loc.computeThreads(0).shellList;
                val shellPrimitives = shellList.getShellPrimitives();
                val numPrimitives = shellList.getNumberOfShellPrimitives();
                val numPairs = shellList.getNumberOfShellPairs();

                var myG : Int = placeId;

                var placePairs:Int = 0;
                var placeInts:Long = 0;
                for(var i:Int=myG; i<numPairs; i++) {
                    if (i == myG) {
                        finish {
                          async myG = G.getAndIncrement();
                          placeInts += comp_loc.computeOneShellPairThreaded(i, numPrimitives, shellPrimitives);
                        }
                        placePairs++;
                    }
                }
                placeTimer.stop(0);
                Console.OUT.printf("\tcompute at %s pairs %i integrals %i %.4g seconds\n", here, placePairs, placeInts, ((placeTimer.total(0) as Double) / 1e9));

                comp_loc.getGMatContributionArray()
            };
            // gather and reduce my gMatrix contribution
            val gatherTimer = new Timer(1);
            gatherTimer.start(0);
            atomic super.cellAdd(placeContribution);
            gatherTimer.stop(0);
            Console.OUT.printf("\tgather from %i %.3g seconds\n", placeId, ((gatherTimer.total(0) as Double) / 1e9));
        }
    }

    /** Compute class for the new code - multi place version */
    static class ComputePlace {
        var density:Density;
        val mol:Molecule[QMAtom];

        val gMatrixContribution:DenseMatrix;
        val jMatrixContribution:DenseMatrix;
        val kMatrixContribution:DenseMatrix;

        private val qCut:DenseMatrix;
        private val dCut:DenseMatrix;
        private val maxEst:Double;
        private var thresh2:Double;

        public val computeThreads = new Array[ComputeThread](Runtime.NTHREADS);

        public def this(N : Int, mol:Molecule[QMAtom], bfs:BasisFunctions, qCut:DenseMatrix, dCut:DenseMatrix, maxEst:Double) {
            this.mol = mol;
            this.qCut = qCut;
            this.dCut = dCut;
            this.maxEst = maxEst;

            gMatrixContribution = new DenseMatrix(N, N);
            jMatrixContribution = new DenseMatrix(N, N);
            kMatrixContribution = new DenseMatrix(N, N);

            // TODO not required with broadcast ateach above - XTENLANG-1725
            val noOfElectrons = mol.getNumberOfElectrons();
            val noOfOccupancies = noOfElectrons / 2;
            density = new Density(N, noOfOccupancies);

            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i) = new ComputeThread(N, bfs, qCut, dCut);
            }
        }

        /**
         * Prepares to calculate G Matrix using the provided
         * density matrix.  Sets partial contribution to G
         * and other temporary arrays to zero.
         */
        public def reset(density : Density) {
            this.density = density;
            gMatrixContribution.reset();
            jMatrixContribution.reset();
            kMatrixContribution.reset();
            recomputeDCut();
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i).reset(density, thresh2);
            }
        }

        /**
         * Prepares to calculate G Matrix using a broadcast
         * of the density matrix from place 0.  Sets partial 
         * contribution to G and other temporary arrays to zero.
         */
        public def reset() {
            // TODO should use ateach broadcast - XTENLANG-1725
            val densityMatrix = density.d;
            Team.WORLD.bcast[Double](here.id, 0, densityMatrix, 0, densityMatrix, 0, densityMatrix.size);
            gMatrixContribution.reset();
            jMatrixContribution.reset();
            kMatrixContribution.reset();
            recomputeDCut();
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
                computeThreads(i).reset(density, thresh2);
            }
        }

        /**
         * Combines the gmat matrix contributions from all places.
         * TODO should use reduce, but not currently implemented in x10.util.Team
         */
        public def allreduceGMat() {
            val gmat = getGMatContributionArray().d;
            Team.WORLD.allreduce[Double](here.id, gmat, 0, gmat, 0, gmat.size, Team.ADD);
        }

        private def recomputeDCut() {
            // Häser & Ahlrichs -- eqn 7  
            val noOfAtoms = mol.getNumberOfAtoms();

            dCut.reset();
            var maxDen:Double = 0.0;
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
                            val bAng  = jbFunc.getMaximumAngularMomentum();
                            val aAng  = iaFunc.getMaximumAngularMomentum();
                            val NA:Int = (aAng+2)*(aAng+1)/2+iaFunc.intIndex;
                            val NB:Int = (bAng+2)*(bAng+1)/2+jbFunc.intIndex;
                            var newMaxDen:Double=0.0;
                            for (iden in iaFunc.intIndex..(NA-1)) for (jden in jbFunc.intIndex..(NB-1))
                                newMaxDen = Math.max(newMaxDen, Math.abs(density(iden,jden)));
                            dCut(iaFunc.intIndex,jbFunc.intIndex) = newMaxDen;
                            if (newMaxDen > maxDen) maxDen = newMaxDen;
                        }
                    }
                 }
             }
            if (Place.FIRST_PLACE == here) {
                Console.OUT.printf("\tmaxDen %.4g\n", maxDen);
            }
            thresh2 = THRESH/maxEst/maxDen;
        }

        /**
         * Computes this place's contribution to the G Matrix
         * given previously calculated J, K matrices.
         */
        public def getGMatContributionArray() {
            computeJMatContribution();
            computeKMatContribution();
             
            for([x,y] in 0..(gMatrixContribution.M-1)*0..(gMatrixContribution.N-1)) {
               gMatrixContribution(x,y) = jMatrixContribution(x,y) - (0.25*kMatrixContribution(x,y));
            }
            return gMatrixContribution;
        }

        public def setValue(a:Int, b:Int, c:Int, d:Int, i:Int, j:Int, k:Int, l:Int) : Boolean {
            val iFunc = mol.getAtom(a).getBasisFunctions().get(i);
            val jFunc = mol.getAtom(b).getBasisFunctions().get(j);
            val kFunc = mol.getAtom(c).getBasisFunctions().get(k);
            val lFunc = mol.getAtom(d).getBasisFunctions().get(l);

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

        public def compute(start:Int, end:Int):Long {
            var totInt:Long = 0;
            val noOfAtoms = mol.getNumberOfAtoms();

            for(var a:Int=start; a<end; a++) {
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

                                       // TODO: 
                                       // Console.OUT.println(a + ", " + b + ", " + c + ", " + d + " | " + i + ", " + j + ", " + k + ", " + l);
                                       totInt += computeThreads(0).computeSingle(iaFunc, jbFunc, kcFunc, ldFunc, density);
                                   } // end l
                               } // centre d
                           } // end k
                       } // centre c
                   } // end j
               } // centre b
             } // end i
          } // centre a
            return totInt;
        }

        public def computeShells(nPairs:Int):Long {
            var totInt:Long = 0;
            val computeThread = computeThreads(0);
            val shellList = computeThread.shellList;
            val shellPrimitives = shellList.getShellPrimitives();
            val numPrimitives = shellList.getNumberOfShellPrimitives();
            val numPairs = shellList.getNumberOfShellPairs();

            for(var i:Int=here.id; i<numPairs; i+=Place.MAX_PLACES) {
                totInt += computeThread.computeOneShellPair(i, numPrimitives, shellPrimitives);
            }
            return totInt;
        }

        private def computeJMatContribution() : DenseMatrix {
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
               jMatrixContribution.cellAdd(computeThreads(i).getJMat());
            }

            return jMatrixContribution;             
        }

        private def computeKMatContribution() : DenseMatrix {
            for(var i:Int=0; i<Runtime.NTHREADS; i++) {
               kMatrixContribution.cellAdd(computeThreads(i).getKMat());
            }

            return kMatrixContribution;
        }

        public def computeOneShellPairThreaded(shellPair:Int, nPrimitives:Int, shellPrimitives:Rail[ContractedGaussian]):Long {
            val a = shellPair / nPrimitives;
            val b = shellPair % nPrimitives;

            val aFunc = shellPrimitives(a);
            val bFunc = shellPrimitives(b);
            return computeOneShellPairThreaded(aFunc, bFunc, nPrimitives, shellPrimitives);
        }

        public def computeOneShellPairThreaded(aFunc:ContractedGaussian, bFunc:ContractedGaussian, nPrimitives:Int, shellPrimitives:Rail[ContractedGaussian]):Long {
            var totInt:Long = 0;

            val aStrt = aFunc.intIndex;
            val bStrt = bFunc.intIndex;
            val aAng  = aFunc.getMaximumAngularMomentum();
            val bAng  = bFunc.getMaximumAngularMomentum();

            val aa = aStrt + aAng;
            val bb = bStrt + bAng;

            if (aa < bb) return 0;

            // Cut-off
            if (qCut(aFunc.intIndex,bFunc.intIndex)<thresh2) return 0;

            val radiusABSquared = aFunc.distanceSquaredFrom(bFunc);

            val currentTriplet = new AtomicInteger(computeThreads.size);
            finish for ([threadNum] in computeThreads) async {
                val computeThread = computeThreads(threadNum);

                var threadInts:Long = 0;

                var myTriplet:Int = threadNum;
                var i:Int = -1;
                for (c in 0..(nPrimitives-1)) {
                    val cFunc = shellPrimitives(c);
                    val cStrt = cFunc.intIndex;
                    val cAng  = cFunc.getMaximumAngularMomentum();
                    val cc = cStrt + cAng;

                    if (aa >= cc) {
                        i++; // a triplet to actually be computed
                        if (i == myTriplet) {
                          threadInts += computeThread.computeOneShellTriplet(aFunc, bFunc, cFunc, radiusABSquared, nPrimitives, shellPrimitives, qCut, dCut);
                          myTriplet = currentTriplet.getAndIncrement();
                        }
                    }
                }
                atomic totInt += threadInts;
            }
            return totInt;
        }
    }

    static class ComputeThread {
        var density:Density;
        var computing:Boolean = false;

        var i:ContractedGaussian, j:ContractedGaussian, 
            k:ContractedGaussian, l:ContractedGaussian;

        val twoEI:TwoElectronIntegrals;
        val shellList:ShellList;
        val jMat:DenseMatrix, kMat:DenseMatrix;
        private val qCut:DenseMatrix;
        private val dCut:DenseMatrix;
        private var thresh2:Double;

        def this(N: Int, bfs:BasisFunctions, qCut:DenseMatrix, dCut:DenseMatrix) {
            val shellList = bfs.getShellList();
            this.twoEI = new TwoElectronIntegrals(shellList.getMaximumAngularMomentum(), bfs.getNormalizationFactors(), THRESH);
            this.shellList = shellList;
            this.qCut = qCut;
            this.dCut = dCut;

            jMat = new DenseMatrix(N,N);
            kMat = new DenseMatrix(N,N);
        }

        public def reset(density:Density, thresh2:Double) {
            this.density = density;
            jMat.reset();
            kMat.reset();
            this.thresh2 = thresh2;
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
                                 density : Density):Long {
            return twoEI.compute2EAndRecord(i, j, k, l, 
                                     shellList,
                                     jMat, kMat, density);
        }

        public def computeSingle2(i:ContractedGaussian, j:ContractedGaussian,
                                  k:ContractedGaussian, l:ContractedGaussian, 
                                  radiusABSquared:Double,
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int,
                                  density : Density):Long {
            return twoEI.compute2EAndRecord2(i, j, k, l,
                                      shellList,
                                      jMat, kMat, density,
                                      radiusABSquared,
                                      aAng, bAng, cAng, dAng, angMomAB,
                                      aStrt, bStrt, cStrt, dStrt,
                                      aLim, bLim);
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

        public def computeOneShellPair(shellPair:Int, nPrimitives:Int, shellPrimitives:Rail[ContractedGaussian]):Long {
            val a = shellPair / nPrimitives;
            val b = shellPair % nPrimitives;

            val aFunc = shellPrimitives(a);
            val bFunc = shellPrimitives(b);
            return computeOneShellPair(aFunc, bFunc, nPrimitives, shellPrimitives);
        }

        public def computeOneShellPair(aFunc:ContractedGaussian, bFunc:ContractedGaussian, nPrimitives:Int, shellPrimitives:Rail[ContractedGaussian]):Long {
            var totInt:Long = 0;

            val aStrt = aFunc.intIndex;
            val bStrt = bFunc.intIndex;
            val aAng  = aFunc.getMaximumAngularMomentum();
            val bAng  = bFunc.getMaximumAngularMomentum();

            val aa = aStrt + aAng;
            val bb = bStrt + bAng;

            if (aa < bb) return 0;

            // Cut-off
            if (qCut(aFunc.intIndex,bFunc.intIndex)<thresh2) return 0;

            val angMomAB = aAng + bAng;
            val aLim = ((aAng+1)*(aAng+2)/2);
            val bLim = ((bAng+1)*(bAng+2)/2);

            val radiusABSquared = aFunc.distanceSquaredFrom(bFunc);

            for (c in 0..(nPrimitives-1)) {
                val cFunc = shellPrimitives(c);
                val cStrt = cFunc.intIndex;
                val cAng  = cFunc.getMaximumAngularMomentum();
                val cc = cStrt + cAng;

                if (aa >= cc) {
                    for (d in 0..(nPrimitives-1)) {
                        val dFunc = shellPrimitives(d);
                        val dStrt = dFunc.intIndex;
                        val dAng  = dFunc.getMaximumAngularMomentum();
                        val dd = dStrt + dAng;

                        if (cc >= dd) {
                            val maxDCut1 = Math.max(dCut(aFunc.intIndex,bFunc.intIndex),dCut(cFunc.intIndex,dFunc.intIndex));
                            val maxDCut2 = .25*Math.max(dCut(cFunc.intIndex,bFunc.intIndex),dCut(aFunc.intIndex,dFunc.intIndex));
                            val maxDCut3 = .25*Math.max(dCut(cFunc.intIndex,aFunc.intIndex),dCut(bFunc.intIndex,dFunc.intIndex));
                            val maxDCut4 = Math.max(maxDCut2,maxDCut3);
                            val maxDCut = Math.max(maxDCut4,maxDCut1);

                            if (qCut(aFunc.intIndex,bFunc.intIndex)*qCut(cFunc.intIndex,dFunc.intIndex)*maxDCut<THRESH)
                                continue;

                            totInt += computeSingle2(aFunc, bFunc, cFunc, dFunc,
                                        radiusABSquared, aAng, bAng, cAng, dAng, angMomAB,
                                        aStrt, bStrt, cStrt, dStrt, aLim, bLim,
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
            return totInt;
        }

        public def computeOneShellTriplet(aFunc:ContractedGaussian, bFunc:ContractedGaussian, cFunc:ContractedGaussian, radiusABSquared:Double, nPrimitives:Int, shellPrimitives:Rail[ContractedGaussian], qCut:DenseMatrix, dCut:DenseMatrix):Long {
            val aStrt = aFunc.intIndex;
            val aAng  = aFunc.getMaximumAngularMomentum();
            val bStrt = bFunc.intIndex;
            val bAng  = bFunc.getMaximumAngularMomentum();

            val cStrt = cFunc.intIndex;
            val cAng  = cFunc.getMaximumAngularMomentum();
            val cc = cStrt + cAng;

            val angMomAB = aAng + bAng;
            val aLim = ((aAng+1)*(aAng+2)/2);
            val bLim = ((bAng+1)*(bAng+2)/2);

            var totInt:Long = 0;
            for (d in 0..(nPrimitives-1)) {
                val dFunc = shellPrimitives(d);
                val dStrt = dFunc.intIndex;
                val dAng  = dFunc.getMaximumAngularMomentum();
                val dd = dStrt + dAng;
                if (cc >= dd) {
                    val maxDCut1 = Math.max(dCut(aStrt,bStrt),dCut(cStrt,dStrt));
                    val maxDCut2 = .25*Math.max(dCut(cStrt,bStrt),dCut(aStrt,dStrt));
                    val maxDCut3 = .25*Math.max(dCut(cStrt,aStrt),dCut(bStrt,dStrt));
                    val maxDCut4 = Math.max(maxDCut2,maxDCut3);
                    val maxDCut = Math.max(maxDCut4,maxDCut1);

                    if (qCut(aFunc.intIndex,bFunc.intIndex)*qCut(cStrt,dStrt)*maxDCut<THRESH)
                        continue;

                    totInt += computeSingle2(aFunc, bFunc, cFunc, dFunc,
                                radiusABSquared, aAng, bAng, cAng, dAng, angMomAB,
                                aStrt, bStrt, cStrt, dStrt, aLim, bLim,
                                density);
                }
            }
            return totInt;
        }

        public def getJMat() = jMat;
        public def getKMat() = kMat;
    }
}

