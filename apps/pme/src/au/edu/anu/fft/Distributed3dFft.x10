/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2011.
 */
package au.edu.anu.fft;

import x10.compiler.NativeCPPInclude;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;

import x10.util.ArrayList;
import x10.util.Team;
import edu.mit.fftw.FFTW;

/**
 * This class implements a distributed three-dimensional FFT using
 * FFTW at the base level for sequential 1D FFTs.
 */
@NativeCPPInclude("FFTW_typedef.h")
public class Distributed3dFft {
    /** The length of one dimension of the 3D source and target arrays. */
    public val dataSize : Int;

    // TODO these arrays should be parameters, but can't be because of GC problems
    private val source : DistArray[Complex](3); 
    private val target : DistArray[Complex](3){self.dist==source.dist};
    private val temp : DistArray[Complex](3){self.dist==source.dist};

    public val dribble = false;

    public def this(dataSize : Int,
                    source : DistArray[Complex](3), 
                    target : DistArray[Complex](3){self.dist==source.dist},
                    temp : DistArray[Complex](3){self.dist==source.dist}) {
        this.dataSize = dataSize;
        this.source = source;
        this.target = target;
        this.temp = temp;
    }

    /**
     * Do a 3D FFT of the source array, storing the result in target.
     * This is implemented as three sequential 1D FFTs, swapping data
     * amongst all places to re-orient the array for each dimension.
     * Makes use of a temp array, which may be the same as either source
     * or target arrays.
     * TODO the third transpose can be avoided IF BdotC is calculated 
     * in the transposed format.  This is not currently the case.
     * This operation would have to be renamed "doFFT3dTranspose"
     * to indicate that the target array has its dimensions transposed.
     */
    public def doFFT3d(forward : Boolean) {
        if (Place.numPlaces()==1L) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target, forward);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            finish ateach(p in Dist.makeUnique(source.dist.places())) {
                // 'scratch' rails, for use in the 1D FFTs
                val oneDSource = new Rail[Complex](dataSize);
                val oneDTarget = new Rail[Complex](dataSize);
                if (dribble) {
                    do1DFftAndTranspose(source, target, oneDSource, oneDTarget, forward);
                    Team.WORLD.barrier();
                    do1DFftAndTranspose(target, temp, oneDSource, oneDTarget, forward);
                    Team.WORLD.barrier();
                    do1DFftAndTranspose(temp, target, oneDSource, oneDTarget, forward);
                } else {
                    val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, oneDSource, oneDTarget, forward);
                    do1DFftToTemp(source, oneDSource, oneDTarget, plan);
                    transposeTempToTarget();
                    Team.WORLD.barrier();
                    do1DFftToTemp(target, oneDSource, oneDTarget, plan);
                    Team.WORLD.barrier();
                    transposeTempToTarget();
                    Team.WORLD.barrier();
                    do1DFftToTemp(target, oneDSource, oneDTarget, plan);
                    Team.WORLD.barrier();
                    transposeTempToTarget();
                    FFTW.fftwDestroyPlan(plan);
                }
            }
/*
 *  Formerly, the above was:
 *
            val oneDSource = DistArray.make[Rail[Complex]](Dist.makeUnique(), (Point) => new Rail[Complex](dataSize));
            val oneDTarget = DistArray.make[Rail[Complex]](oneDSource.dist, (Point) => new Rail[Complex](dataSize));
            finish ateach(p1 in oneDSource) do1DFftToTemp(source, oneDSource(p1), oneDTarget(p1), forward);
            finish ateach(p1 in oneDSource) transposeTempToTarget();
            finish ateach(p1 in oneDSource) do1DFftToTemp(target, oneDSource(p1), oneDTarget(p1), forward);
            finish ateach(p1 in oneDSource) transposeTempToTarget();
            finish ateach(p1 in oneDSource) do1DFftToTemp(target, oneDSource(p1), oneDTarget(p1), forward);
            finish ateach(p1 in oneDSource) transposeTempToTarget();

 * or replace with ScalableTreeBarrier once XTENLANG-1660 is resolved
            val c = Clock.make();
            finish {
                for (p1 in source.dist.places()) async(p1) clocked (c) {
                    // 'scratch' rails, for use in the 1D FFTs
                    val oneDSource = Rail.make[Complex](dataSize);
                    val
                    do1DFftToTemp(source, oneDSource, oneDTarget, forward);
                    transposeTempToTarget();
                    next;
                    do1DFftToTemp(target, oneDSource, oneDTarget, forward);
                    next;
                    transposeTempToTarget();
                    next;
                    do1DFftToTemp(target, oneDSource, oneDTarget, forward);
                    next;
                    transposeTempToTarget();
                }
                c.drop();
            }
*/
        }
    }

    /**
     * Performs a 1D FFT for each 1D slice along the Z dimension,
     * and store the result in the temp array.
     */
    private def do1DFftToTemp(source : DistArray[Complex](3),
                                     oneDSource : Rail[Complex],
                                     oneDTarget : Rail[Complex],
                                     plan : FFTW.FFTWPlan) {
        val localSource = source.getLocalPortion();
        val localRegion = localSource.region as Region(3){rect};
        val localTemp = temp.getLocalPortion();
        val gridRegionWithoutZ = (localRegion.eliminate(2)) as Region(2){rect};
        for ([i,j] in gridRegionWithoutZ) {
            // TODO need to copy into Array - can use raw()?
            for (k in 0..(dataSize-1)) {
                oneDSource(k) = localSource(i,j,k);
            }
            FFTW.fftwExecute(plan);
            for (k in 0..(dataSize-1)) {
                localTemp(i,j,k) = oneDTarget(k);
            }
        }
    }

    /**
     * Shuffles array around all places by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays contain complete
     * pencils in the second dimension.
     */
    private def transposeTempToTarget() {
        val sourceDist = temp.dist | here;
        val sourceStartX = sourceDist.region.min(0);
        val sourceEndX = sourceDist.region.max(0);
        val sourceStartY = sourceDist.region.min(1);
        val sourceEndY = sourceDist.region.max(1);
        val localTemp = temp.getLocalPortion();
        val target = this.target; // TODO shouldn't be necessary XTENLANG-1913;
        finish {
            for (p2 in temp.dist.places()) {
                val targetDist = temp.dist | p2;
                val startX = Math.max(sourceStartX, targetDist.region.min(1));
                val endX = Math.min(sourceEndX, targetDist.region.max(1));
                val startY = Math.max(sourceStartY, targetDist.region.min(2));
                val endY = Math.min(sourceEndY, targetDist.region.max(2));
                val startZ = targetDist.region.min(0);
                val endZ = targetDist.region.max(0);

                val transferRegion = Region.makeRectangular([(startX..endX), (startY..endY), (startZ..endZ)]);
                if (transferRegion.size() > 0) {
                    if (p2 == here) {
                        val localTarget = target.getLocalPortion();
                        // do a synchronous in-place transpose from temp to target
                        for ([x,y,z] in transferRegion) {
                            localTarget(z,x,y) = localTemp(x,y,z);
                        }
                    } else {
                        // collect elements to transfer to remote place
                        val elementsToTransfer = new Rail[Complex](transferRegion.size());
                        var i : Int = 0;
                        for ([x,y,z] in transferRegion) {
                            elementsToTransfer(i++) = localTemp(x,y,z);
                        }
                        at(p2) async {
                            val localTarget = target.getLocalPortion();
                            var i2 : Int = 0;
                            for ([x,y,z] in transferRegion) {
                                // transpose dimensions
                                localTarget(z,x,y) = elementsToTransfer(i2++);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Performs a 1D FFT for each 1D pencil along the Z dimension,
     * and dribbles the results to be transposed in the target array.
     * Assumes block,block distributed data, that is, each place
     * holds a rectangular block of data over some portion of the zeroth
     * and first dimensions, that is complete in the second dimension.
     */
    private def do1DFftAndTranspose(source : DistArray[Complex](3),
                                    target : DistArray[Complex](3){self.dist==source.dist},
                                    oneDSource : Rail[Complex],
                                    oneDTarget : Rail[Complex],
                                    forward : Boolean) {
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, oneDSource, oneDTarget, forward);
        val myRegion = source.dist(here);
        val gridRegionWithoutZ = (myRegion.eliminate(2)) as Region(2){rect};
        finish for ([i,j] in gridRegionWithoutZ) {
            // TODO need to copy into Rail - can use raw()?
            for (k in 0..(dataSize-1)) {
                oneDSource(k) = source(i,j,k);
            }
            FFTW.fftwExecute(plan);
            // "Dribble" this pencil to appropriate target places
            for (p2 in target.dist.places()) {
                val targetRegion = target.dist(p2);
                if (i >= targetRegion.min(1) && i <= targetRegion.max(1)) {
                    val startX = targetRegion.min(0);
                    val endX = targetRegion.max(0);
                    val numToTransfer = endX-startX+1;
                    if (numToTransfer > 0) {
                        val elementsToTransfer = new Rail[Complex](numToTransfer);
                        var k:Int = 0;
                        for (x in startX..endX) {
                            elementsToTransfer(k++) = oneDTarget(x);
                        }
                        at(p2) async {
                            var x:Int = 0;
                            for (k2 in startX..endX) {
                                target(k2,i,j) = elementsToTransfer(x++);
                            }
                        }
                    }
                }
            }
        }
        FFTW.fftwDestroyPlan(plan);
    }
}

