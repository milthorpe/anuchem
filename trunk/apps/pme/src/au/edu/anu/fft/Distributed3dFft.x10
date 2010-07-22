package au.edu.anu.fft;

import x10.util.GrowableRail;
import edu.mit.fftw.FFTW;

/**
 * This class implements a distributed three-dimensional FFT using
 * FFTW at the base level for sequential 1D FFTs.
 */
public class Distributed3dFft {
    /** The length of one dimension of the 3D source and target arrays. */
    public global val dataSize : Int;

    // TODO these arrays should be parameters, but can't be because of GC problems
    private global val source : DistArray[Complex](3); 
    private global val target : DistArray[Complex](3){self.dist==source.dist};
    private global val temp : DistArray[Complex](3){self.dist==source.dist};

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
    public global def doFFT3d(forward : Boolean) {
        if (source.dist.constant) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target, forward);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            do1DFftToTemp(source, forward);
            transposeTempToTarget();
            do1DFftToTemp(target, forward);
            transposeTempToTarget();
            do1DFftToTemp(target, forward);
            transposeTempToTarget();
        }
    }

    /**
     * Performs a 1D FFT for each 1D slice along the Z dimension,
     * and store the result in the temp array.
     */
    private global def do1DFftToTemp(source : DistArray[Complex](3),
                                     forward : Boolean) {
        finish ateach ((p1) in Dist.makeUnique(source.dist.places())) {
            val oneDSource = Rail.make[Complex](dataSize);
            val oneDTarget = Rail.make[Complex](dataSize);
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, oneDSource, oneDTarget, forward);
            val mySource = source.dist | here;
            val gridRegionWithoutZ = (mySource.region().eliminate(2)) as Region(2);
            for ((i,j) in gridRegionWithoutZ) {
                // TODO need to copy into ValRail - can use raw()?
                for(var k : Int = 0; k < dataSize; k++) {
                    oneDSource(k) = source(i,j,k);
                }
                FFTW.fftwExecute(plan);
                for(var k : Int = 0; k < dataSize; k++) {
                    temp(i,j,k) = oneDTarget(k);
                }
            }
            FFTW.fftwDestroyPlan(plan);
        }
    }

    /**
     * Shuffles array around all places by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays contain complete
     * slabs or lines in the second dimension.
     */
    private global def transposeTempToTarget() {
        finish ateach (p1 in Dist.makeUnique(temp.dist.places())) {
            val sourceDist = temp.dist | here;
            val sourceStartX = sourceDist.region.min(0);
            val sourceEndX = sourceDist.region.max(0);
            val sourceStartY = sourceDist.region.min(1);
            val sourceEndY = sourceDist.region.max(1);
            foreach (p2 in temp.dist.places()) {
                val targetDist = temp.dist | p2;
                val startX = Math.max(sourceStartX, targetDist.region.min(1));
                val endX = Math.min(sourceEndX, targetDist.region.max(1));
                val startY = Math.max(sourceStartY, targetDist.region.min(2));
                val endY = Math.min(sourceEndY, targetDist.region.max(2));
                val startZ = targetDist.region.min(0);
                val endZ = targetDist.region.max(0);

                val sourcePoints : Region(3) = [startX..endX, startY..endY, startZ..endZ];
                val elementsToTransfer = Rail.make[Complex](sourcePoints.size());
                var i : Int = 0;
                for (p in sourcePoints) {
                    elementsToTransfer(i++) = temp(p);
                }
                val toTransfer = elementsToTransfer as ValRail[Complex];
                at (p2) {
                    val targetPoints : Region(3) = [startX..endX, startY..endY, startZ..endZ];
                    var i : Int = 0;
                    for ((x,y,z) in targetPoints) {
                        target(z,x,y) = toTransfer(i++);
                    }
                }
            }
        }
    }
}
