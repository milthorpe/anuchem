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
        if (source.dist.constant) {
            // all source data at one place.  use local 3D FFT rather than distributed
            val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(dataSize, dataSize, dataSize, source, target, forward);
            FFTW.fftwExecute(plan);
            FFTW.fftwDestroyPlan(plan); 
        } else {
            doFFTForOneDimension(source, temp, forward);
            transposeArray(temp, target);
            doFFTForOneDimension(target, temp, forward);
            transposeArray(temp, target);
            doFFTForOneDimension(target, temp, forward);
            transposeArray(temp, target);
        }
    }

    /**
     * Performs a 1D FFT for each 1D slice along the Z dimension.
     */
    private global def doFFTForOneDimension(source : DistArray[Complex](3), 
                                     target : DistArray[Complex](3){self.dist==source.dist},
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
                    target(i,j,k) = oneDTarget(k);
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
    private global def transposeArray(source : DistArray[Complex](3), 
                             target : DistArray[Complex](3){self.dist==source.dist}) {
        finish ateach (p1 in Dist.makeUnique(source.dist.places())) {
            val sourceDist = source.dist | here;
            val sourceStart = sourceDist.region.min(0);
            val sourceEnd = sourceDist.region.max(0);
            foreach (p2 in source.dist.places()) {
                val targetDist = source.dist | p2;
                val targetStart = targetDist.region.min(0);
                val targetEnd = targetDist.region.max(0);
                
                val toTransfer = ValRail.make[Complex]((targetEnd - targetStart + 1) * (sourceEnd - sourceStart + 1) * dataSize, (n : Int) => source(mapPoint(n,sourceStart,targetStart,targetEnd)));
                at (p2) {transpose(toTransfer, target, sourceStart, sourceEnd, targetStart, targetEnd);}
            }   
        }
    }

    private global safe def mapPoint(n : Int, iStart : Int, jStart : Int, jEnd : Int) : Point(3) {
        val JSIZE = jEnd - jStart + 1;
        val ISIZE = dataSize * JSIZE;
        val i = n / (ISIZE);
        val j = (n - i * ISIZE) / JSIZE;
        val k = n - i * ISIZE - j * JSIZE;
        return Point.make(i+iStart,j,k+jStart);
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
        val KSIZE = targetEnd - targetStart + 1;
        val JSIZE = dataSize * KSIZE;
        for ((j) in sourceStart..sourceEnd) {
            val jOffset = (j-sourceStart) * JSIZE;
            for ((k) in 0..dataSize-1) {
                val kOffset = k * KSIZE;
                for ((i) in targetStart..targetEnd) {
                    val offset = kOffset + jOffset + (i-targetStart);
                    //Console.OUT.println(Point.make(i,j,k) + " <= " + offset);
                    target(i,j,k) = source(offset);
                }
            }
        }
    }
}
