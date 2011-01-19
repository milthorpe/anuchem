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

import x10.compiler.Native;
import x10.util.ArrayList;
import x10.util.Team;
import edu.mit.fftw.FFTW;

/**
 * This class implements a distributed three-dimensional FFT using
 * FFTW at the base level for sequential 1D FFTs.
 */
public class Distributed3dFft {
    /** The length of one dimension of the 3D source and target arrays. */
    public val dataSize : Int;

    // TODO these arrays should be parameters, but can't be because of GC problems
    private val source : DistArray[Complex](3); 
    private val target : DistArray[Complex](3){self.dist==source.dist};
    private val temp : DistArray[Complex](3){self.dist==source.dist};

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
        if (Place.MAX_PLACES==1) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target, forward);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            finish {
                for (p1 in source.dist.places()) async at(p1) {
                    // 'scratch' rails, for use in the 1D FFTs
                    val oneDSource = new Array[Complex](dataSize);
                    val oneDTarget = new Array[Complex](dataSize);
                    do1DFftToTemp(source, oneDSource, oneDTarget, forward);
                    transposeTempToTarget();
                    Team.WORLD.barrier(here.id);
                    do1DFftToTemp(target, oneDSource, oneDTarget, forward);
                    Team.WORLD.barrier(here.id);
                    transposeTempToTarget();
                    Team.WORLD.barrier(here.id);
                    do1DFftToTemp(target, oneDSource, oneDTarget, forward);
                    Team.WORLD.barrier(here.id);
                    transposeTempToTarget();
                }
            }
/*
 *  Formerly, the above was:
 *
            val oneDSource = DistArray.make[Array[Complex](1){rail}](Dist.makeUnique(), (Point) => new Array[Complex](dataSize));
            val oneDTarget = DistArray.make[Array[Complex](1){rail}](oneDSource.dist, (Point) => new Array[Complex](dataSize));
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
                                     oneDSource : Array[Complex](1){rail},
                                     oneDTarget : Array[Complex](1){rail},
                                     forward : Boolean) {
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, oneDSource, oneDTarget, forward);
        val mySource = source.dist | here;
        val gridRegionWithoutZ = (mySource.region().eliminate(2)) as Region(2){rect};
        for ([i,j] in gridRegionWithoutZ) {
            // TODO need to copy into Array - can use raw()?
            for ([k] in 0..(dataSize-1)) {
                oneDSource(k) = source(i,j,k);
            }
            FFTW.fftwExecute(plan);
            for ([k] in 0..(dataSize-1)) {
                temp(i,j,k) = oneDTarget(k);
            }
        }
        FFTW.fftwDestroyPlan(plan);
    }

    /**
     * Shuffles array around all places by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays contain complete
     * slabs or lines in the second dimension.
     */
    private def transposeTempToTarget() {
        val sourceDist = temp.dist | here;
        val sourceStartX = sourceDist.region.min(0);
        val sourceEndX = sourceDist.region.max(0);
        val sourceStartY = sourceDist.region.min(1);
        val sourceEndY = sourceDist.region.max(1);
        finish {
            for (p2 in temp.dist.places()) {
                val targetDist = temp.dist | p2;
                val startX = Math.max(sourceStartX, targetDist.region.min(1));
                val endX = Math.min(sourceEndX, targetDist.region.max(1));
                val startY = Math.max(sourceStartY, targetDist.region.min(2));
                val endY = Math.min(sourceEndY, targetDist.region.max(2));
                val startZ = targetDist.region.min(0);
                val endZ = targetDist.region.max(0);

                val transferRegion = (startX..endX) * (startY..endY) * (startZ..endZ);
                if (transferRegion.size() > 0) {
                    val elementsToTransfer = new Array[Complex](transferRegion.size());
                    var i : Int = 0;
                    for ([x,y,z] in transferRegion) {
                        elementsToTransfer(i++) = temp(x,y,z);
                    }
                    async at (p2) {
                        var i : Int = 0;
                        for ([x,y,z] in transferRegion) {
                            // transpose dimensions
                            target(z,x,y) = elementsToTransfer(i++);
                        }
                    }
                }
            }
        }
    }
}

