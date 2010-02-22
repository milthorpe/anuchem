package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;
import edu.mit.fftw.FFTW;
import x10.array.BaseArray;

public class PME {
    /** The number of grid lines in each dimension of the simulation unit cell. */
    public val gridSize : ValRail[Int](3);

    /** Double representations of the various grid dimensions */
    public val K1 : Double;
    public val K2 : Double; 
    public val K3 : Double;

    /** The edges of the unit cell. */
    public val edges : ValRail[Vector3d](3);

    /** The conjugate reciprocal vectors for each dimension. */
    public val edgeReciprocals : ValRail[Vector3d](3);

    global val gridRegion : Region(3);

	val atoms : ValRail[MMAtom];
    
    /** The order of B-spline interpolation */
    val splineOrder : Int;

    /** The Ewald coefficient beta */
    val beta : Double;

    /** The direct sum cutoff distance in Angstroms */
    val cutoff : Double;

    /** Translation vectors for neighbouring unit cells (the 26 cells surrounding the origin cell) */
    val imageTranslations : Array[Vector3d](3);

    // TODO should be shared local to calculateEnergy()
    var directEnergy : Double = 0.0;
    var directSum : Double = 0.0;
    var selfEnergy: Double = 0.0;
    var correctionEnergy : Double = 0.0; // TODO masklist

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
            atoms : ValRail[MMAtom],
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
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        val imageTranslationRegion = [-1..1,-1..1,-1..1] as Region(3);
        this.imageTranslations = Array.make[Vector3d](imageTranslationRegion, (p(i,j,k) : Point(3)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k)));
        Console.OUT.println("PME for " + atoms.length + " particles.");
        Console.OUT.println("Box edges: " + edges + " volume: " + getVolume());
        Console.OUT.println("Grid size: " + gridSize);
        Console.OUT.println("spline order: " + splineOrder + " Beta: " + beta + " Cutoff: " + cutoff);
    }
	
    public def calculateEnergy() : Double {
        // TODO do we need multi-threaded FFTW, or multiple activities (with FFT) per place?
        FFTW.fftwInitThreads();
        FFTW.fftwPlanWithNThreads(Runtime.INIT_THREADS);

        finish foreach ((i) in 0..atoms.length-1) {
            var myDirectEnergy : Double = 0.0;
            var myDirectSum : Double = 0.0;
            // NOTE include i==j as this contributes image components
            for (var j : Int = 0; j < atoms.length; j++) {
                val rjri = Vector3d(atoms(j).centre.sub(atoms(i).centre as Tuple3d));
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

        Console.OUT.println("directSum / selfEnergy");
        
        val Q = getGriddedCharges();
        Console.OUT.println("Q");
        /*
        for (p in Q) {
            if (Q(p) != Complex.ZERO) {
                Console.OUT.println(p + " = " + Q(p));
            }
        }
        */

        //val Qinv = DFT.dft3D(Q, false);

        val Qinv = Array.make[Complex](gridRegion);
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan3d(gridSize(0), gridSize(1), gridSize(2), Q as BaseArray[Complex], Qinv as BaseArray[Complex], false);
        FFTW.fftwExecute(plan);
        FFTW.fftwDestroyPlan(plan);
        
        Console.OUT.println("Qinv");

        val B = getBArray();
        Console.OUT.println("B");
        val C = getCArray();
        Console.OUT.println("C");

        val BdotC = Array.make[Double](gridRegion, (p : Point(gridRegion.rank)) => B(p) * C(p));
        Console.OUT.println("BdotC");

        val thetaRecConvQInv = Array.make[Complex](gridRegion, (p : Point(gridRegion.rank)) => BdotC(p) * Qinv(p));
        val thetaRecConvQ = Array.make[Complex](gridRegion);

        val plan2 : FFTW.FFTWPlan = FFTW.fftwPlan3d(gridSize(0), gridSize(1), gridSize(2), thetaRecConvQInv as BaseArray[Complex], thetaRecConvQ as BaseArray[Complex], true);
        FFTW.fftwExecute(plan2);
        FFTW.fftwDestroyPlan(plan2);
        
        //val thetaRecConvQ = DFT.dft3D(thetaRecConvQInv, true);
        /*for (p in thetaRecConvQ) {
            if (thetaRecConvQ(p) != Complex.ZERO) {
                Console.OUT.println(p + " = " + thetaRecConvQ(p));
            }
        }*/
        Console.OUT.println("thetaRecConvQ");

        val reciprocalEnergy = getReciprocalEnergy(Q, thetaRecConvQ);

        Console.OUT.println("directEnergy = " + directEnergy);
        Console.OUT.println("directSum = " + directSum);
        Console.OUT.println("selfEnergy = " + selfEnergy);
        Console.OUT.println("correctionEnergy = " + correctionEnergy);
        Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val total = directSum + reciprocalEnergy.re + (correctionEnergy + selfEnergy);
        val error = directEnergy - total;
        Console.OUT.println("error = " + error + " relative error = " + Math.abs(error) / Math.abs(total));
        FFTW.fftwCleanupThreads();
        return directSum + reciprocalEnergy.re + (correctionEnergy + selfEnergy);
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     */
    public def getGriddedCharges() : Array[Complex](3){self.region==gridRegion} {
        val Q = Array.make[Double](gridRegion);
        finish foreach ((i) in 0..atoms.length-1) {
            val atom = atoms(i);
            val q = atom.charge;
            val u = getScaledFractionalCoordinates(Vector3d(atom.centre as Tuple3d));
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

                        atomic {
                            Q((k1 + gridSize(0)) % gridSize(0),
                              (k2 + gridSize(1)) % gridSize(1),
                              (k3 + gridSize(2)) % gridSize(2)) += gridPointContribution;
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
        val Qcomplex = Array.make[Complex](gridRegion, (m : Point(gridRegion.rank)) => Complex(Q(m), 0.0));
        return Qcomplex;
    }

    /**
     * @param Q the gridded charge array as defined in Eq. 4.6
     * @param thetaRecConvQ the reciprocal pair potential as defined in eq. 4.7
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(Q : Array[Complex]{self.region==gridRegion}, thetaRecConvQ : Array[Complex]{self.region==gridRegion}) {
        var reciprocalEnergy : Complex = Complex.ZERO;
        for (m in gridRegion) {
            val gridPointContribution = Q(m) * thetaRecConvQ(m);
            if (gridPointContribution.re != 0.0) {
                //Console.OUT.println("gridPointContribution( " + m + ") = " + gridPointContribution);
            }
            reciprocalEnergy = reciprocalEnergy + gridPointContribution;
        }
        reciprocalEnergy = reciprocalEnergy / 2.0;
        return reciprocalEnergy;
    }

    /**
     * @return the array B as defined by Eq 4.8 and 4.4
     */
    public def getBArray() {
        val B = Array.make[Double](gridRegion);
        // TODO for (m(m1,m2,m3) in gridRegion) {
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
                    var sumK3 : Complex = Complex.ZERO;
                    for ((k) in 0..(splineOrder-2)) {
                        sumK3 = sumK3 + bSpline4(k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
                    }
                    val b3 = (Math.exp(2.0 * Math.PI * (splineOrder - 1.0) * m3D / K3 * Complex.I) / sumK3).abs();
                    //Console.OUT.println("b1 = " + b1 + " b2 = " + b2 + " b3 = " + b3);
                    B(m1,m2,m3) = b1 * b1 * b2 * b2 * b3 * b3;
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
        val C = Array.make[Double](gridRegion);
        val V = getVolume();
        //Console.OUT.println("V = " + V);
        finish foreach (m(m1,m2,m3) in gridRegion) {
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
        C(0,0,0) = 0.0;
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
    public def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
        return new Vector3d(edgeReciprocals(0).mul(gridSize(0)).dot(r), edgeReciprocals(1).mul(gridSize(1)).dot(r), edgeReciprocals(2).mul(gridSize(2)).dot(r));
    }

    /**
     * Gets the volume V of the unit cell.
     */
    public safe def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }
}
