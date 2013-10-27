/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2013.
 */
package au.edu.anu.fft;

import x10.compiler.NativeCPPInclude;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;

import x10.util.Team;
import edu.mit.fftw.FFTW;

/**
 * This class implements a distributed three-dimensional FFT using
 * FFTW at the base level for sequential 1D FFTs.
 */
@NativeCPPInclude("FFTW_typedef.h")
public class DistributedReal3dFft_SPMD {
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
    public static def doFFT3d(source:DistArray[Double](3), target:DistArray[Complex](3){self.dist==source.dist}, temp:DistArray[Complex](3){self.dist==source.dist}) {
        val dataSize = (source.region.max(0)-source.region.min(0)) as Int + 1n;
        if (Place.MAX_PLACES==1L) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan);
        } else {
            val localSource = source.getLocalPortion();
            val localRegion = localSource.region as Region(3){rect};
            val gridRegionWithoutZ = (localRegion.eliminate(2)) as Region(2){rect};
            val howMany = gridRegionWithoutZ.size() as Int;

            val planR2C:FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, howMany, source, temp);
            FFTW.fftwExecute(planR2C);
            // FFTW only computes half-complex DFT.  fill in the other half
            val diagonal = dataSize/2n+1n;
            for ([i,j] in gridRegionWithoutZ) {
                for (k in diagonal..(dataSize-1n)) {
                    temp(i,j,k) = temp(i,j,dataSize-k).conjugate();
                }
            }
            transpose[Complex](temp, target);
            Team.WORLD.barrier();
            val planC2C:FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, howMany, target, temp, true);
            FFTW.fftwExecute(planC2C);
            Team.WORLD.barrier();
            transpose[Complex](temp, target);
            Team.WORLD.barrier();
            FFTW.fftwExecute(planC2C);
            Team.WORLD.barrier();
            transpose[Complex](temp, target);
            FFTW.fftwDestroyPlan(planR2C);
            FFTW.fftwDestroyPlan(planC2C);
        }
    }

    public static def doFFT3d(source:DistArray[Complex](3), target:DistArray[Double](3){self.dist==source.dist}, temp:DistArray[Complex](3){self.dist==source.dist}, temp2:DistArray[Double](3){self.dist==source.dist}) {
        val dataSize = (source.region.max(0)-source.region.min(0)) as Int + 1n;
        if (Place.MAX_PLACES==1L) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            val localTarget = source.getLocalPortion();
            val localRegion = localTarget.region as Region(3){rect};
            val gridRegionWithoutZ = (localRegion.eliminate(2)) as Region(2){rect};
            val howMany = gridRegionWithoutZ.size() as Int;

            val planC2C:FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, howMany, source, temp, false);
            FFTW.fftwExecute(planC2C);
            Team.WORLD.barrier();
            transpose[Complex](temp, source);
            Team.WORLD.barrier();
            FFTW.fftwExecute(planC2C);
            Team.WORLD.barrier();
            transpose[Complex](temp, source);
            Team.WORLD.barrier();
            val planC2R:FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, howMany, source, temp2);
            FFTW.fftwExecute(planC2R);
            transpose[Double](temp2, target);
            FFTW.fftwDestroyPlan(planC2C);
            FFTW.fftwDestroyPlan(planC2R);
        }
    }

    /**
     * Shuffles array around all places by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays contain complete
     * pencils in the second dimension.
     */
    private static def transpose[T](source:DistArray[T](3), target:DistArray[T](3){self.dist==target.dist}){T haszero} {
        val sourceDist = source.dist | here;
        val sourceStartX = sourceDist.region.min(0);
        val sourceEndX = sourceDist.region.max(0);
        val sourceStartY = sourceDist.region.min(1);
        val sourceEndY = sourceDist.region.max(1);
        val localSource = source.getLocalPortion();
        finish {
            // transpose elements to every place, starting with next place
            // and ending with 'here'
            var p2:Place = here.next();
            do {
                val targetDist = source.dist | p2;
                val startX = Math.max(sourceStartX, targetDist.region.min(1));
                val endX = Math.min(sourceEndX, targetDist.region.max(1));
                val startY = Math.max(sourceStartY, targetDist.region.min(2));
                val endY = Math.min(sourceEndY, targetDist.region.max(2));
                val startZ = targetDist.region.min(0);
                val endZ = targetDist.region.max(0);

                val transferRegion = Region.make(startX..endX, startY..endY, startZ..endZ);
                if (transferRegion.size() > 0) {
                    if (p2 == here) {
                        val localTarget = target.getLocalPortion();
                        // synchronous in-place transpose from source to target
                        for ([x,y,z] in transferRegion) {
                            localTarget(z,x,y) = localSource(x,y,z);
                        }
                    } else {
                        // collect elements to transfer to remote place
                        val elementsToTransfer = new Rail[T](transferRegion.size());
                        var i : Long = 0;
                        for ([x,y,z] in transferRegion) {
                            elementsToTransfer(i++) = localSource(x,y,z);
                        }
                        at(p2) async {
                            val localTarget = target.getLocalPortion();
                            var i2:Long = 0;
                            for ([x,y,z] in transferRegion) {
                                // transpose dimensions
                                localTarget(z,x,y) = elementsToTransfer(i2++);
                            }
                        }
                    }
                }
                p2 = p2.next();
            } while (p2 != here.next());
        }
    }
}

