package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;
import edu.mit.fftw.FFTW;
import au.edu.anu.util.Timer;

public class PME {
    // TODO enum - XTENLANG-1118
    public const TIMER_INDEX_TOTAL : Int = 0;
    public const TIMER_INDEX_DIRECT : Int = 1;
    public const TIMER_INDEX_GRIDCHARGES : Int = 2;
    public const TIMER_INDEX_INVFFT : Int = 3;
    public const TIMER_INDEX_THETARECCONVQ : Int = 4;
    public const TIMER_INDEX_RECIPROCAL : Int = 5;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(6);

    /** The number of grid lines in each dimension of the simulation unit cell. */
    private global val gridSize : ValRail[Int](3);

    /** Double representations of the various grid dimensions */
    private global val K1 : Double;
    private global val K2 : Double; 
    private global val K3 : Double;

    /** The edges of the unit cell. */
    private global val edges : ValRail[Vector3d](3);

    /** The conjugate reciprocal vectors for each dimension. */
    private global val edgeReciprocals : ValRail[Vector3d](3);

    private global val gridRegion : Region(3);
    private global val gridDist : Dist(3);
    
    /** The order of B-spline interpolation */
    private global val splineOrder : Int;

    /** The Ewald coefficient beta */
    private global val beta : Double;

    /** The direct sum cutoff distance in Angstroms */
    private global val cutoff : Double;

    /** Translation vectors for neighbouring unit cells (the 26 cells surrounding the origin cell) */
    private global val imageTranslations : Array[Vector3d](3);

	private val atoms : ValRail[MMAtom!];

    private val B : Array[Double]{self.dist==gridDist};
    private val C : Array[Double]{self.dist==gridDist};
    private val BdotC : Array[Double]{self.dist==gridDist};

    // TODO should be shared local to calculateEnergy()
    private var directEnergy : Double = 0.0;
    private var directSum : Double = 0.0;
    private var selfEnergy: Double = 0.0;
    private var correctionEnergy : Double = 0.0; // TODO masklist
    private var reciprocalEnergy : Double = 0.0;

    /**
     * Creates a new particle mesh Ewald method.
     * @param edges the edge vectors of the unit cell
     * @param gridSize the number of grid lines in each dimension of the unit cell
     * @param atoms the atoms in the unit cell
     * @param splineOrder the order n of B-spline interpolation
     * @param beta the Ewald coefficient beta
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(edges : ValRail[Vector3d](3),
            gridSize : ValRail[Int](3),
            atoms : ValRail[MMAtom!],
            splineOrder : Int,
            beta : Double,
            cutoff : Double) {
        this.gridSize = gridSize;
        K1 = gridSize(0) as Double;
        K2 = gridSize(1) as Double;
        K3 = gridSize(2) as Double;
        this.edges = edges;
        this.edgeReciprocals = ValRail.make[Vector3d](3, (i : Int) => edges(i).inverse());
        this.atoms = atoms;
        val r = Region.makeRectangular(0, gridSize(0)-1);
        gridRegion = (r * [0..(gridSize(1)-1)] * [0..(gridSize(2)-1)]) as Region(3);
        gridDist = Dist.makeBlock(gridRegion, 0);
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        val imageTranslationRegion = [-1..1,-1..1,-1..1] as Region(3);
        this.imageTranslations = Array.make[Vector3d](imageTranslationRegion, (p(i,j,k) : Point(3)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k)));

        Console.OUT.println("PME for " + atoms.length + " particles.");
        Console.OUT.println("Box edges: " + edges + " volume: " + getVolume());
        Console.OUT.println("Grid size: " + gridSize);
        Console.OUT.println("spline order: " + splineOrder + " Beta: " + beta + " Cutoff: " + cutoff);
        Console.OUT.println("gridDist = " + gridDist);

        B = getBArray();
        C = getCArray();
        BdotC = Array.make[Double](gridDist, (p : Point(gridDist.region.rank)) => B(p) * C(p));
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        timer.start(TIMER_INDEX_DIRECT);
        finish foreach ((i) in 0..atoms.length-1) {
            var myDirectEnergy : Double = 0.0;
            var myDirectSum : Double = 0.0;
            // NOTE include i==j as this contributes image components
            for (var j : Int = 0; j < atoms.length; j++) {
                val rjri = new Vector3d(atoms(j).centre.sub(atoms(i).centre as Tuple3d));
                // rough (non-Euclidean, 1D) distance cutoff to avoid unnecessary distance calculations
                // for (p(n1,n2,n3) in imageTranslations) {
                for (var n1:Int = -1; n1<=1; n1++) {
                    for (var n2:Int = -1; n2<=1; n2++) {
                        for (var n3:Int = -1; n3<=1; n3++) {
                            if (! (i==j && (n1 | n2 | n3) == 0)) {
                                val imageDistance = rjri.add(imageTranslations(n1,n2,n3)).length();
                                val chargeProduct = atoms(i).charge * atoms(j).charge;
                                myDirectEnergy += chargeProduct / imageDistance;
                                if (imageDistance < cutoff) {
                                    val imageDirectComponent = chargeProduct * Math.erfc(beta * imageDistance) / imageDistance;
                                    //Console.OUT.println("imageDistance = " + imageDistance + " imageDirectComponent = " + imageDirectComponent);
                                    myDirectSum += imageDirectComponent;
                                    //Console.OUT.println("distance = " + distance + " directEnergy component = " + chargeProduct / distance);
                                    
                                }
                            }
                        }
                    }
                }
            }
            atomic {
                directEnergy += myDirectEnergy;
                directSum += myDirectSum;
                selfEnergy += atoms(i).charge * atoms(i).charge;
            }
        }
        selfEnergy = -beta / Math.sqrt(Math.PI) * selfEnergy;
        directEnergy = directEnergy / 2.0;
        directSum = directSum / 2.0;
        timer.stop(TIMER_INDEX_DIRECT);
        
        val Q = getGriddedCharges();

        timer.start(TIMER_INDEX_INVFFT);
        val temp = Array.make[Complex](gridDist); 
        val Qinv = Array.make[Complex](gridDist);
        doFFT3d(Q, Qinv, temp, false);
        timer.stop(TIMER_INDEX_INVFFT);

        timer.start(TIMER_INDEX_THETARECCONVQ);
        val thetaRecConvQInv = Array.make[Complex](gridDist, (p : Point(gridDist.region.rank)) => BdotC(p) * Qinv(p));
        val thetaRecConvQ = Array.make[Complex](gridDist);
        doFFT3d(thetaRecConvQInv, thetaRecConvQ, temp, true);
        timer.stop(TIMER_INDEX_THETARECCONVQ);

        reciprocalEnergy = getReciprocalEnergy(Q, thetaRecConvQ);

        //Console.OUT.println("directEnergy = " + directEnergy);
        //Console.OUT.println("directSum = " + directSum);
        //Console.OUT.println("selfEnergy = " + selfEnergy);
        //Console.OUT.println("correctionEnergy = " + correctionEnergy);
        //Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val totalEnergy = directSum + reciprocalEnergy + (correctionEnergy + selfEnergy);
        val error = directEnergy - totalEnergy;
        Console.OUT.println("error = " + error + " relative error = " + Math.abs(error) / Math.abs(totalEnergy));

        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     */
    public def getGriddedCharges() : Array[Complex](3){self.dist==gridDist} {
        timer.start(TIMER_INDEX_GRIDCHARGES);
        val Q = Array.make[Double](gridDist);
        finish foreach ((i) in 0..atoms.length-1) {
            val atom = atoms(i);
            val q = atom.charge;
            val u = getScaledFractionalCoordinates(new Vector3d(atom.centre as Tuple3d));
            //Console.OUT.println("atom( " + i + " ) charge = " + q + " coords = " + u);
            //var atomContribution : Double = 0.0;

            val u1c = Math.ceil(u.i) as Int;
            val u2c = Math.ceil(u.j) as Int;
            val u3c = Math.ceil(u.k) as Int;
            for (var k1 : Int = u1c - splineOrder; k1 < u1c; k1++) {
                for (var k2 : Int = u2c - splineOrder; k2 < u2c; k2++) {
                    for (var k3 : Int = u3c - splineOrder; k3 < u3c; k3++) {
                        val gridPointContribution = q
                                     * bSpline4(u.i - k1)
                                     * bSpline4(u.j - k2)
                                     * bSpline4(u.k - k3);

                        val p = Point.make((k1 + gridSize(0)) % gridSize(0),
                                           (k2 + gridSize(1)) % gridSize(1),
                                           (k3 + gridSize(2)) % gridSize(2));
                        async(Q.dist(p)) {
                            atomic {
                                Q(p) += gridPointContribution;
                            }
                        }
                        /*
                        Console.OUT.println("Q(" + (k1 + gridSize(0)) % gridSize(0) + "," + 
                                                   (k2 + gridSize(1)) % gridSize(1) + "," + 
                                                   (k3 + gridSize(2)) % gridSize(2) + 
                        ") += " + gridPointContribution);
                        atomContribution += gridPointContribution;
                        */
                    }
                }
            }
            //Console.OUT.println("atomContribution = " + atomContribution);
        }
        val Qcomplex = Array.make[Complex](gridDist, (m : Point(gridDist.region.rank)) => Complex(Q(m), 0.0));
        timer.stop(TIMER_INDEX_GRIDCHARGES);
        return Qcomplex;
    }

    /**
     * @param Q the gridded charge array as defined in Eq. 4.6
     * @param thetaRecConvQ the reciprocal pair potential as defined in eq. 4.7
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(Q : Array[Complex]{self.dist==gridDist}, thetaRecConvQ : Array[Complex]{self.dist==gridDist}) {
        timer.start(TIMER_INDEX_RECIPROCAL);
        finish for ((p1) in Dist.makeUnique(Place.places)) {
            async (Place.places(p1)) {
                var myReciprocalEnergy : Double = 0.0;
                // TODO single-place parallel
                for ((i,j,k) in gridDist | here) {
                    val gridPointContribution = Q(i,j,k) * thetaRecConvQ(i,j,k);
                    myReciprocalEnergy = myReciprocalEnergy + gridPointContribution.re;
                }
                myReciprocalEnergy = myReciprocalEnergy / 2.0;
                val myReciprocalEnergyFinal = myReciprocalEnergy;
                at (this) {atomic{reciprocalEnergy += myReciprocalEnergyFinal;}};
            }
        }
        timer.stop(TIMER_INDEX_RECIPROCAL);
        return reciprocalEnergy;
    }

    /**
     * @return the array B as defined by Eq 4.8 and 4.4
     */
    public def getBArray() {
        val B = Array.make[Double](gridDist);
        // TODO for (m(m1,m2,m3) in gridDist) {
        finish foreach ((m1) in 0..gridSize(0)-1) {
            val m1D = m1 as Double;
            var sumK1 : Complex = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK1 = sumK1 + bSpline4(k+1) * Math.exp(2.0 * Math.PI * m1D * k / K1 * Complex.I);
            }
            val b1 = (Math.exp(2.0 * Math.PI * (splineOrder - 1.0) * m1D / K1* Complex.I) / sumK1).abs();

            for (var m2 : Int = 0; m2 < gridSize(1); m2++) {
                val m2D = m2 as Double;
                var sumK2 : Complex = Complex.ZERO;
                for ((k) in 0..(splineOrder-2)) {
                    sumK2 = sumK2 + bSpline4(k+1) * Math.exp(2.0 * Math.PI * m2D * k / K2 * Complex.I);
                }
                val b2 = (Math.exp(2.0 * Math.PI * (splineOrder - 1.0) * m2D / K2 * Complex.I) / sumK2).abs();
                
                for (var m3 : Int = 0; m3 < gridSize(2); m3++) {
                    // TODO smarter distribution
                    val m3D = m3 as Double;
                    val m = Point.make(m1,m2,m3);
                    at(B.dist(m)) {
                        var sumK3 : Complex = Complex.ZERO;
                        for ((k) in 0..(splineOrder-2)) {
                            sumK3 = sumK3 + bSpline4(k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
                        }
                        val b3 = (Math.exp(2.0 * Math.PI * (splineOrder - 1.0) * m3D / K3 * Complex.I) / sumK3).abs();
                        //Console.OUT.println("b1 = " + b1 + " b2 = " + b2 + " b3 = " + b3);
                        B(m) = b1 * b1 * b2 * b2 * b3 * b3;
                    }
                    //Console.OUT.println("B(" + m1 + "," + m2 + "," + m3 + ") = " + B(m1,m2,m3));
                }
            }
        }
        return B;
    }

    /**
     * @return the array C as defined by Eq 3.9
     */
    public def getCArray() {
        val C = Array.make[Double](gridDist);
        val V = getVolume();
        //Console.OUT.println("V = " + V);
        finish ateach (m(m1,m2,m3) in gridDist) {
            val m1prime = m1 <= K1/2 ? m1 : m1 - K1;
            val m2prime = m2 <= K2/2 ? m2 : m2 - K2;
            val m3prime = m3 <= K3/2 ? m3 : m3 - K3;
            val mVec = edgeReciprocals(0).mul(m1prime).add(edgeReciprocals(1).mul(m2prime)).add(edgeReciprocals(2).mul(m3prime));
            //Console.OUT.println("mVec = " + mVec);
            val mSquared = mVec.dot(mVec);
            //Console.OUT.println("mSquared = " + mSquared);
            //Console.OUT.println("numerator (exp) = " + Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)));
            //Console.OUT.println("denominator = " + (mSquared * Math.PI * V));
            C(m) = Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)) / (mSquared * Math.PI * V);
            //Console.OUT.println("C" + m + " = " + C(m));
        }
        at (C.dist(0,0,0)) {
            C(0,0,0) = 0.0;
        }
        return C;
    }

    /* 
     * Gets the nth order B-spline M_n(u) as per Eq. 4.1
     */
    public static safe def bSpline(n : Int, u : Double) : Double {
        if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 2) {
            return 1.0 - Math.abs(u - 1.0);
        } else {
            return u / (n - 1) * bSpline(n-1, u) + (n - u) / (n - 1) * bSpline(n-1, u-1.0);
        }
    }

    /* 
     * Gets the 4th order B-spline M_4(u) as per Eq. 4.1
     */
    private static safe def bSpline4(u : Double) : Double {
        if (u <= 0.0 || u >= 4) {
            return 0.0;
        } else {
            return u / 3 * bSpline3(u) + (4 - u) / 3 * bSpline3(u-1.0);
        }
    }

    /* 
     * Gets the 3rd order B-spline M_3(u) as per Eq. 4.1
     */
    private static safe def bSpline3(u : Double) : Double {
        if (u <= 0.0 || u >= 3) {
            return 0.0;
        } else {
            return u / 2 * bSpline2(u) + (3 - u) / 2 * bSpline2(u-1.0);
        }
    }

    /* 
     * Gets the 2nd order B-spline M_2(u) as per Eq. 4.1
     */
    private static safe def bSpline2(u : Double) : Double {
        if (u <= 0.0 || u >= 2) {
            return 0.0;
        } else {
            return 1.0 - Math.abs(u - 1.0);
        }
    }
    
    /** Gets scaled fractional coordinate u as per Eq. 3.1 */
    public global safe def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
        return new Vector3d(edgeReciprocals(0).mul(gridSize(0)).dot(r), edgeReciprocals(1).mul(gridSize(1)).dot(r), edgeReciprocals(2).mul(gridSize(2)).dot(r));
    }

    /**
     * Gets the volume V of the unit cell.
     */
    private global safe def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }

    /**
     * Do a 3D FFT of the source array, storing the result in target.
     * This is implemented as three sequential 1D FFTs, swapping data
     * amongst all places twice to re-orient the array for each dimension.
     * Makes use of a temp array, which may be the same as either source
     * or target arrays.
     */
    private def doFFT3d(source : Array[Complex](3){self.dist==gridDist}, 
                        target : Array[Complex](3){self.dist==gridDist},
                        temp : Array[Complex](3){self.dist==gridDist},
                        forward : Boolean) {
        doFFTForOneDimension(source, temp, forward);
        if (forward) {
            shuffleArray(temp, target);
        } else {
            shuffleArrayReverse(temp, target);
        }
        doFFTForOneDimension(target, temp, forward);
        if (forward) {
            shuffleArray(temp, target);
        } else {
            shuffleArrayReverse(temp, target);
        }
        doFFTForOneDimension(target, target, forward);
    }

    /**
     * Does a 1D FFT for each 1D slice along the first dimension.
     */
    private def doFFTForOneDimension(source : Array[Complex](3){self.dist==gridDist}, 
                                     target : Array[Complex](3){self.dist==gridDist},
                                     forward : Boolean) {
        finish for ((p1) in Dist.makeUnique(Place.places)) {
            async (Place.places(p1)) {
                val size = gridSize(1);
                val oneDSource = Rail.make[Complex](size);
                val oneDTarget = Rail.make[Complex](size);
                val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(size, oneDSource, oneDTarget, forward);
                val mySource = gridDist | here;
                val gridRegionWithoutFirst = (mySource.region().projection(0) * mySource.region().projection(2)) as Region(2);
                for ((i,k) in gridRegionWithoutFirst) {
                    // TODO need to copy into ValRail?
                    for(var j : Int = 0; j < size; j++) {
                        oneDSource(j) = source(i,j,k);
                    }
                    FFTW.fftwExecute(plan);
                    for(var j : Int = 0; j < size; j++) {
                        target(i,j,k) = oneDTarget(j);
                    }
                }
                FFTW.fftwDestroyPlan(plan);
            }
        }
    }

    /**
     * "Shuffles" array around all places by transposing the zeroth dimension
     * to the second, the second to the first and the first to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays are block 
     * distributed along the zeroth dimension.
     */
    public def shuffleArray(source : Array[Complex](3){self.dist==gridDist}, 
                             target : Array[Complex](3){self.dist==gridDist}) {
        finish for ((p1) in Dist.makeUnique(Place.places)) {
            async (Place.places(p1)) {
                val mySource = gridDist | here;
                foreach ((p2) in Dist.makeUnique(Place.places)) {
                    val place2ZeroDimension = (gridDist | Place.places(p2)).region.projection(0);
                    val myContribution = mySource.region().projection(0) * place2ZeroDimension * mySource.region.projection(2) as Region(3);
                    for ((i,j,k) in myContribution) {
                        val s = source(i,j,k);
                        at(Place.places(p2)) { target(j,k,i) = s;};
                    }
                }
            }
        }
    }

    /**
     * "Shuffles" array around all places "in reverse" by transposing the zeroth dimension
     * to the first, the first to the second and the second to the zeroth.
     * Assumes NxNxN arrays, and that source and target arrays are block 
     * distributed along the zeroth dimension.
     */
    public def shuffleArrayReverse(source : Array[Complex](3){self.dist==gridDist}, 
                             target : Array[Complex](3){self.dist==gridDist}) {
        for ((p1) in Dist.makeUnique(Place.places)) {
            at (Place.places(p1)) {
                val mySource = gridDist | here;
                for ((p2) in Dist.makeUnique(Place.places)) {
                    val place2ZeroDimension = (gridDist | Place.places(p2)).region.projection(0);
                    val myContribution = mySource.region().projection(0) * mySource.region.projection(1) * place2ZeroDimension  as Region(3);
                    for ((i,j,k) in myContribution) {
                        val s = source(i,j,k);
                        at(Place.places(p2)) { target(k,i,j) = s;};
                    }
                }
            }
        }
    }
}
