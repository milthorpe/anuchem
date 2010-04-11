package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.fft.Distributed3dFft;
import au.edu.anu.util.Timer;

/**
 * This class implements a Smooth Particle Mesh Ewald method to calculate
 * the potential of a system of charged particles, as described in
 * Essmann et al. "A Smooth Particle Mesh Ewald method", J. Comp. Phys. 101,
 * pp.8577-8593 (1995) DOI: 10.1063/1.470117
 */
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
    private global val imageTranslations : Array[Vector3d](3){rect};

	private val atoms : ValRail[MMAtom!];

    private global val fft : Distributed3dFft;

    private global val B : DistArray[Double]{self.dist==gridDist};
    private global val C : DistArray[Double]{self.dist==gridDist};
    private global val BdotC : DistArray[Double]{self.dist==gridDist};

    // TODO should be shared local to calculateEnergy() - XTENLANG-404
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
        fft = new Distributed3dFft(gridSize(0));
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
        val imageTranslationRegion = [-1..1,-1..1,-1..1] as Region(3){rect};
        this.imageTranslations = new Array[Vector3d](imageTranslationRegion, (p(i,j,k) : Point(3)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k)));

        Console.OUT.println("PME for " + atoms.length + " particles.");
        Console.OUT.println("Box edges: " + edges + " volume: " + getVolume());
        Console.OUT.println("Grid size: " + gridSize);
        Console.OUT.println("spline order: " + splineOrder + " Beta: " + beta + " Cutoff: " + cutoff);
        Console.OUT.println("gridDist = " + gridDist);

        B = getBArray();
        C = getCArray();
        BdotC = DistArray.make[Double](gridDist);
        // TODO Array.lift not implemented XTENLANG-376
        // BdotC = B.lift((b:Double,c:Double)=>b*c, C);
	    finish ateach(p in gridDist) {
		    BdotC(p) = B(p) * C(p);
	    }
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        timer.start(TIMER_INDEX_DIRECT);
        finish foreach ((i) in 0..atoms.length-1) {
            // TODO assume cutoff < edgeLength, so we don't need i==j
            var myDirectSum : Double = 0.0;
            for (var j : Int = 0; j < i; j++) {
                val rjri = new Vector3d(atoms(j).centre.sub(atoms(i).centre as Tuple3d));
                // TODO rough (non-Euclidean, 1D) distance cutoff to avoid unnecessary distance calculations
                for ((n1,n2,n3) in imageTranslations) {
                    if (! (i==j && (n1 | n2 | n3) == 0)) {
                        val imageDistance = rjri.add(imageTranslations(n1,n2,n3)).length();
                        if (imageDistance < cutoff) {
                            val chargeProduct = atoms(i).charge * atoms(j).charge;
                            val imageDirectComponent = chargeProduct * Math.erfc(beta * imageDistance) / imageDistance;
                            myDirectSum += imageDirectComponent;                                    
                        }
                    }
                }
            }
            
            // TODO this is slow because of lack of optimized atomic - XTENLANG-321
            atomic {
                directSum += myDirectSum * 2;
                selfEnergy += atoms(i).charge * atoms(i).charge;
            }
        }
        selfEnergy = -beta / Math.sqrt(Math.PI) * selfEnergy;
        directSum = directSum / 2.0;
        timer.stop(TIMER_INDEX_DIRECT);

        //Console.OUT.println("directSum");
        
        val Q = getGriddedCharges();

        timer.start(TIMER_INDEX_INVFFT);
        val temp = DistArray.make[Complex](gridDist); 
        val Qinv = DistArray.make[Complex](gridDist);
        fft.doFFT3d(Q, Qinv, temp, false);
        timer.stop(TIMER_INDEX_INVFFT);

        //Console.OUT.println("Qinv");

        timer.start(TIMER_INDEX_THETARECCONVQ);
        val thetaRecConvQInv = DistArray.make[Complex](gridDist, (p : Point(gridDist.region.rank)) => BdotC(p) * Qinv(p));
        val thetaRecConvQ = DistArray.make[Complex](gridDist);
        fft.doFFT3d(thetaRecConvQInv, thetaRecConvQ, temp, true);
        timer.stop(TIMER_INDEX_THETARECCONVQ);

        //Console.OUT.println("thetaRecConvQ");

        reciprocalEnergy = getReciprocalEnergy(Q, thetaRecConvQ);

        //Console.OUT.println("directEnergy = " + directEnergy);
        //Console.OUT.println("directSum = " + directSum);
        //Console.OUT.println("selfEnergy = " + selfEnergy);
        //Console.OUT.println("correctionEnergy = " + correctionEnergy);
        //Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val totalEnergy = directSum + reciprocalEnergy + (correctionEnergy + selfEnergy);

        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     */
    public def getGriddedCharges() : DistArray[Complex](3){self.dist==gridDist} {
        timer.start(TIMER_INDEX_GRIDCHARGES);
        val Q = DistArray.make[Double](gridDist);
        finish foreach ((i) in 0..atoms.length-1) {
            val atom = atoms(i);
            val q = atom.charge;
            val u = getScaledFractionalCoordinates(new Vector3d(atom.centre as Tuple3d));
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
                            // TODO this is slow because of lack of optimized atomic - XTENLANG-321
                            atomic { Q(p) += gridPointContribution; }
                        }
                    }
                }
            }
        }
        val Qcomplex = DistArray.make[Complex](gridDist, (m : Point(gridDist.region.rank)) => Complex(Q(m), 0.0));
        timer.stop(TIMER_INDEX_GRIDCHARGES);
        return Qcomplex;
    }

    /**
     * @param Q the gridded charge array as defined in Eq. 4.6
     * @param thetaRecConvQ the reciprocal pair potential as defined in eq. 4.7
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(Q : DistArray[Complex]{self.dist==gridDist}, thetaRecConvQ : DistArray[Complex]{self.dist==gridDist}) {
        timer.start(TIMER_INDEX_RECIPROCAL);
        finish ateach ((p1) in Dist.makeUnique(gridDist.places())) {
            var myReciprocalEnergy : Double = 0.0;
            // TODO single-place parallel - requires efficient atomic XTENLANG-321
            for ((i,j,k) in gridDist | here) {
                val gridPointContribution = Q(i,j,k) * thetaRecConvQ(i,j,k);
                myReciprocalEnergy = myReciprocalEnergy + gridPointContribution.re;
            }
            myReciprocalEnergy = myReciprocalEnergy / 2.0;
            val myReciprocalEnergyFinal = myReciprocalEnergy;
            at (this) {atomic{reciprocalEnergy += myReciprocalEnergyFinal;}};
        }
        timer.stop(TIMER_INDEX_RECIPROCAL);
        return reciprocalEnergy;
    }

    /**
     * @return the array B as defined by Eq 4.8 and 4.4
     */
    public def getBArray() {
        val B = DistArray.make[Double](gridDist);
        // TODO ateach (m(m1,m2,m3) in gridDist) {
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
                    val m3D = m3 as Double;
                    val m = Point.make(m1,m2,m3);
                    async(B.dist(m)) {
                        var sumK3 : Complex = Complex.ZERO;
                        for ((k) in 0..(splineOrder-2)) {
                            sumK3 = sumK3 + bSpline4(k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
                        }
                        val b3 = (Math.exp(2.0 * Math.PI * (splineOrder - 1.0) * m3D / K3 * Complex.I) / sumK3).abs();
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
        val C = DistArray.make[Double](gridDist);
        val V = getVolume();
        //Console.OUT.println("V = " + V);
        finish ateach (m(m1,m2,m3) in gridDist) {
            val m1prime = m1 <= K1/2 ? m1 : m1 - K1;
            val m2prime = m2 <= K2/2 ? m2 : m2 - K2;
            val m3prime = m3 <= K3/2 ? m3 : m3 - K3;
            val mVec = edgeReciprocals(0).mul(m1prime).add(edgeReciprocals(1).mul(m2prime)).add(edgeReciprocals(2).mul(m3prime));
            val mSquared = mVec.dot(mVec);
            C(m) = Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)) / (mSquared * Math.PI * V);
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
}
