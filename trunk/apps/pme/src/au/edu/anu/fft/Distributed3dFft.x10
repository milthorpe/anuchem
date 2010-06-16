package au.edu.anu.fft;

import edu.mit.fftw.FFTW;

/**
 * This class implements a distributed three-dimensional FFT using
 * FFTW at the base level for sequential 1D FFTs.
 */
public class Distributed3dFft {
    /** The length of one dimension of the 3D source and target arrays. */
    public global val dataSize : Int;

    public def this(dataSize : Int) {
        this.dataSize = dataSize;
        //FFTW.fftwInitThreads();
        //FFTW.fftwPlanWithNThreads(Runtime.INIT_THREADS);
    }

    // TODO a neater way to cleanup the FFTW environment
    public def cleanup() {
        //FFTW.fftwCleanupThreads();
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
    public global def doFFT3d(source : DistArray[Complex](3), 
                        target : DistArray[Complex](3){self.dist==source.dist},
                        temp : DistArray[Complex](3){self.dist==source.dist},
                        forward : Boolean) {
        if (source.dist.onePlace == here) {
            // all source data is here.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target, forward);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            doFFTForOneDimension(source, temp, forward);
            if (forward) {
                transposeArray(temp, target);
            } else {
                transposeArrayReverse(temp, target);
            }
            doFFTForOneDimension(target, temp, forward);
            if (forward) {
                transposeArray(temp, target);
            } else {
                transposeArrayReverse(temp, target);
            }
            doFFTForOneDimension(target, temp, forward);
            if (forward) {
                transposeArray(temp, target);
            } else {
                transposeArrayReverse(temp, target);
            }
        }
    }

    /**
     * Performs a 1D FFT for each 1D slice along the first dimension.
     */
    private global def doFFTForOneDimension(source : DistArray[Complex](3), 
                                     target : DistArray[Complex](3){self.dist==source.dist},
                                     forward : Boolean) {
        finish ateach ((p1) in Dist.makeUnique(source.dist.places())) {
            val oneDSource = Rail.make[Complex](dataSize);
            val oneDTarget = Rail.make[Complex](dataSize);
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(dataSize, oneDSource, oneDTarget, forward);
            val mySource = source.dist | here;
            val gridRegionWithoutFirst = (mySource.region().projection(0) * mySource.region().projection(2)) as Region(2);
            for ((i,k) in gridRegionWithoutFirst) {
                // TODO need to copy into ValRail - can use raw()?
                for(var j : Int = 0; j < dataSize; j++) {
                    oneDSource(j) = source(i,j,k);
                }
                FFTW.fftwExecute(plan);
                for(var j : Int = 0; j < dataSize; j++) {
                    target(i,j,k) = oneDTarget(j);
                }
            }
            FFTW.fftwDestroyPlan(plan);
        }
    }

    /**
     * Shuffles array around all places by transposing the zeroth dimension
     * to the second, the second to the first and the first to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays are block 
     * distributed along the zeroth dimension.
     */
    private global def transposeArray(source : DistArray[Complex](3), 
                             target : DistArray[Complex](3){self.dist==source.dist}) {
        finish ateach ((p1) in Dist.makeUnique(source.dist.places())) {
            val sourceDist = source.dist | here;
            val sourceStart = sourceDist.region.min(0);
            val sourceEnd = sourceDist.region.max(0);
            foreach (p2 in source.dist.places()) {
                val targetDist = source.dist | p2;
                val targetStart = targetDist.region.min(0);
                val targetEnd = targetDist.region.max(0);
                val toTransfer = ValRail.make[Complex]((targetEnd - targetStart + 1) * (sourceEnd - sourceStart + 1) * dataSize, (n : Int) => source(mapPoint(n,sourceStart,sourceEnd,targetStart)));
                at (p2) {transpose(toTransfer, target, sourceStart, sourceEnd, targetStart, targetEnd);}
            }
        }
    }

    /**
     * Shuffles array around all places "in reverse" by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays are block 
     * distributed along the zeroth dimension.
     */
    private global def transposeArrayReverse(source : DistArray[Complex](3), 
                             target : DistArray[Complex](3){self.dist==source.dist}) {
        finish ateach ((p1) in Dist.makeUnique(source.dist.places())) {
            val sourceDist = source.dist | here;
            val sourceStart = sourceDist.region.min(0);
            val sourceEnd = sourceDist.region.max(0);
            foreach (p2 in source.dist.places()) {
                val targetDist = source.dist | p2;
                val targetStart = targetDist.region.min(0);
                val targetEnd = targetDist.region.max(0);
                val toTransfer = ValRail.make[Complex](dataSize * (targetEnd - targetStart + 1) * (sourceEnd - sourceStart + 1), (n : Int) => source(mapPointReverse(n,sourceStart,sourceEnd,targetStart)));
                at (p2) {transposeReverse(toTransfer, target, sourceStart, sourceEnd, targetStart, targetEnd);}
            }
        }
    }

    private global safe def mapPoint(n : Int, iStart : Int, iEnd : Int, jStart : int) : Point(3) {
        val JOFFSET = (iEnd - iStart+1) * dataSize;
        val IOFFSET = dataSize;
        val j = n / JOFFSET;
        val i = (n - j * JOFFSET) / IOFFSET;
        val k = n - j * JOFFSET - i * IOFFSET;
        return Point.make(i+iStart,j+jStart,k);
    }

    private global safe def mapPointReverse(n : Int, iStart : Int, iEnd : Int, kStart : int) : Point(3) {
        val KOFFSET = (iEnd - iStart+1) * dataSize;
        val IOFFSET = dataSize;
        val k = n / KOFFSET;
        val i = (n - k * KOFFSET) / IOFFSET;
        val j = n - k * KOFFSET - i * IOFFSET;
        return Point.make(i+iStart,j,k+kStart);
    }

    /**
     * Copies a chunk of a source array into the target array.
     * The elements are shoehorned into a ValRail of length K * (endI - startI) .
     */
    private global safe def transpose(source : ValRail[Complex],
                             target : DistArray[Complex](3),
                             sourceStart : Int,
                             sourceEnd : Int,
                             targetStart : Int,
                             targetEnd : Int) : Void {
        for ((i) in targetStart..targetEnd) {
            val iOffset = (i-targetStart) * (sourceEnd-sourceStart+1) * dataSize;
            for ((k) in sourceStart..sourceEnd) {
                val kOffset = (k - sourceStart) * dataSize;
                for (var j : Int = 0; j < dataSize; j++) {
                    val offset = kOffset + iOffset + j;
                    target(i,j,k) = source(offset);
                }
            }
        }
    }

    /**
     * Copies a chunk of a source array into the target array.
     * The elements are shoehorned into a ValRail of length K * (endI - startI) .
     */
    private global safe def transposeReverse(source : ValRail[Complex],
                             target : DistArray[Complex](3),
                             sourceStart : Int,
                             sourceEnd : Int,
                             targetStart : Int,
                             targetEnd : Int) : Void {
        for ((i) in targetStart..targetEnd) {
            val iOffset = (i-targetStart) * (sourceEnd-sourceStart+1) * dataSize;
            for ((j) in sourceStart..sourceEnd) {
                val jOffset = (j-sourceStart) * dataSize;
                for (var k : Int = 0; k < dataSize; k++) {
                    val offset = jOffset + iOffset + k;
                    target(i,j,k) = source(offset);
                }
            }
        }
    }
}
