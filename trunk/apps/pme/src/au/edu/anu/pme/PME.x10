package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.fft.DFT;
import au.edu.anu.chem.mm.MMAtom;

public class PME {
    /** The number of grid lines in each dimension of the simulation unit cell. */
    public val gridSize : ValRail[Int](3);

    /** Double representations of the various grid dimensions */
    public val K1 : Double;
    public val K2 : Double; 
    public val K3 : Double;

    /** The edges of the unit cell. TODO: non-cubic cells */
    public val edges : ValRail[Vector3d](3);

    /** The conjugate reciprocal vectors for each dimension. */
    public val edgeReciprocals : ValRail[Vector3d](3);

	val atoms : ValRail[MMAtom];

    val gridRegion : Region(3);
    
    /** The order of B-spline interpolation */
    val splineOrder : Int;

    /** The Ewald coefficient beta */
    val beta : Double;

    /** The direct sum cutoff distance in Angstroms */
    val cutoff : Double;

    /**
     * Creates a new particle mesh Ewald method.
     * @param size the side length of the unit cube // TODO non-cubic boxes
     * @param n the order of B-spline interpolation
     * @param beta the Ewald coefficient beta
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
    }
	
    public def calculateEnergy() : Double {
        Console.OUT.println("PME for " + atoms.length + " particles.");
        Console.OUT.println("Box edges: " + edges + " volume: " + getVolume());
        Console.OUT.println("Grid size: " + gridSize);
        Console.OUT.println("spline order: " + splineOrder + " Beta: " + beta + " Cutoff: " + cutoff);
        var directEnergy : Double = 0.0;
        var directSum : Double = 0.0;
        var selfEnergy: Double = 0.0;
        var correctionEnergy : Double = 0.0;
        for ((i) in 0..(atoms.length - 1)) {
            selfEnergy -= beta / Math.sqrt(Math.PI) * atoms(i).charge * atoms(i).charge;
            for ((j) in 0..(i - 1)) {
                val pairEnergy : Double = atoms(j).pairEnergy(atoms(i));
                directEnergy += 2 * pairEnergy;
                val rjri = new Vector3d(atoms(j).centre.sub(atoms(i).centre));
                val distance = rjri.length();
                Console.OUT.println("beta * distance = " + (beta * distance) + " erf = " + Math.erf(beta * distance));
                correctionEnergy -= atoms(i).charge * atoms(j).charge * Math.erf(beta * distance) / distance;
                for (var n1:Int = -1; n1<=1; n1++) {
                    for (var n2:Int = -1; n2<=1; n2++) {
                        for (var n3:Int = -1; n3<=1; n3++) {
                            if ((n1 & n2 & n3) != 0) {
                                val imageDistance = rjri.add(edges(0).mul(n1)).add(edges(1).mul(n2)).add(edges(2).mul(n3)).length();
                                Console.OUT.println("imageDistance = " + imageDistance);
                                if (imageDistance < cutoff) {
                                    val imageDirectComponent = atoms(i).charge * atoms(j).charge * Math.erfc(beta * imageDistance) / imageDistance;
                                    Console.OUT.println("imageDirectComponent = " + imageDirectComponent);
                                    directSum += imageDirectComponent;
                                }
                            }
                        }
                    }
                }
            }
        }
        Console.OUT.println("directEnergy = " + directEnergy);
        Console.OUT.println("directSum = " + directSum);
        Console.OUT.println("selfEnergy = " + selfEnergy);
        Console.OUT.println("correctionEnergy = " + correctionEnergy);

        val Q = getGriddedCharges();
        Console.OUT.println("Q");
        val Qinv = DFT.dft3D(Q, false);
        Console.OUT.println("Qinv");
        val B = getBArray();
        Console.OUT.println("B");
        val C = getCArray();
        Console.OUT.println("C");
        val BdotC = Array.make[Double](gridRegion, (p : Point(gridRegion.rank)) => B(p) * C(p));
        Console.OUT.println("BdotC");
        val size = (K1*K2*K3);
        val thetaRecConvQ = DFT.dft3D(Array.make[Complex](gridRegion, (p : Point(gridRegion.rank)) => BdotC(p) * Qinv(p)), true);
        Console.OUT.println("thetaRecConvQ");
        /*for (p in thetaRecConvQ) {
            if (thetaRecConvQ(p) != Complex.ZERO) {
                Console.OUT.println(p + " = " + thetaRecConvQ(p));
            }
        }*/
        var reciprocalEnergy : Complex = Complex.ZERO;
        for (m in gridRegion) {
            val gridPointContribution = Q(m) * thetaRecConvQ(m);
            if (gridPointContribution.re != 0.0) {
                Console.OUT.println("gridPointContribution( " + m + ") = " + gridPointContribution);
            }
            reciprocalEnergy = reciprocalEnergy + gridPointContribution;
        }
        reciprocalEnergy = reciprocalEnergy / 2.0;
        Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        return directSum + reciprocalEnergy.re + (correctionEnergy + selfEnergy);
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     */
    public def getGriddedCharges() : Array[Complex](3)!{self.rect,self.zeroBased} {
        val Q = Array.make[Complex](gridRegion, (Point) => Complex.ZERO);
        for ((i) in 0..atoms.length-1) {
            val atom = atoms(i);
            val q = atom.charge;
            val u = getScaledFractionalCoordinates(new Vector3d(atom.centre));
            val u1i = u.i as Int;
            val u2i = u.j as Int;
            val u3i = u.k as Int;
            // TODO: not currently possible to do "wraparound" Region in X10
            for (var kk1:Int=u1i-splineOrder;kk1<=u1i+splineOrder;kk1++) {
                val k1=(kk1+gridSize(0))%gridSize(0);
                for (var kk2:Int=u2i-splineOrder;kk2<=u2i+splineOrder;kk2++) {
                    val k2=(kk2+gridSize(1))%gridSize(1);
                    for (var kk3:Int=u3i-splineOrder;kk3<=u3i+splineOrder;kk3++) {
                        val k3=(kk3+gridSize(2))%gridSize(2);
                //for ((n1, n2, n3) : Point(3) in [-1..1,-1..1,-1..1]) { // TODO: for n > gridSize
                        Q(k1,k2,k3) = Q(k1,k2,k3) + q
                                     * bSpline(splineOrder, (u.i - kk1) % gridSize(0))
                                     * bSpline(splineOrder, (u.j - kk2) % gridSize(1))
                                     * bSpline(splineOrder, (u.k - kk3) % gridSize(2));
                    }
                }
            }
        }
        /*
        // Works exactly as per definition, but horribly inefficient.
        for (k(k1, k2, k3) in gridRegion) {
            var sum : Double = 0.0;
            for ((i) in 0..atoms.length-1) {
                val atom = atoms(i);
                val q = atom.charge;
                val u = getScaledFractionalCoordinates(new Vector3d(atom.centre));
                val u1i = u.i;
                val u2i = u.j;
                val u3i = u.k;
                for ((n1, n2, n3) : Point(3) in [-n..n,-n..n,-n..n]) {
                    sum += q * bSpline(n, u1i - k1 - n1 * gridSize(0))
                             * bSpline(n, u2i - k2 - n2 * gridSize(1))
                             * bSpline(n, u3i - k3 - n3 * gridSize(2));
                 }
            }
            Q(k) = su
        }*/
        return Q;
    }

    /**
     * Calculates the reciprocal pair potential theta_rec used in Eq. 4.7
     */
    public def getReciprocalPairPotential() {

    } 

    /**
     * Calculates the array B as defined by Eq 4.8 and 4.4
     */
    public def getBArray() {
        val B = Array.make[Double](gridRegion);
        for (m(m1,m2,m3) in gridRegion) {
            // TODO: use proper array regions!
            val m1D = m1 as Double;
            val m2D = m2 as Double;
            val m3D = m3 as Double;
            var sumK : Complex = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK = sumK + bSpline(splineOrder, k+1) * (2.0 * Math.PI * m1D * k / K1 * Complex.I).exp();
            }
            val b1 = ((2.0 * Math.PI * (splineOrder - 1.0) * m1D / K1* Complex.I).exp() / sumK).abs();
            sumK = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK = sumK + bSpline(splineOrder, k+1) * (2.0 * Math.PI * m2D * k / K2 * Complex.I).exp();
            }
            val b2 = ((2.0 * Math.PI * (splineOrder - 1.0) * m2D / K2 * Complex.I).exp() / sumK).abs();
            sumK = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK = sumK + bSpline(splineOrder, k+1) * (2.0 * Math.PI * m3D * k / K3 * Complex.I).exp();
            }
            val b3 = ((2.0 * Math.PI * (splineOrder - 1.0) * m3D / K3 * Complex.I).exp() / sumK).abs();
            B(m) = b1 * b1 * b2 * b2 * b3 * b3;
        }
        return B;
    }

    /**
     * Calculates the array C as defined by Eq 3.9
     */
    public def getCArray() {
        val C = Array.make[Double](gridRegion);
        val V = getVolume();
        Console.OUT.println("V = " + V);
        for (m(m1,m2,m3) in gridRegion) {
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
    public def bSpline(n : Int, u : Double) : Double {
        if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 2) {
            return 1.0 - Math.abs(u - 1.0);
        } else {
            return u / (n - 1) * bSpline(n-1, u) + (n - u) / (n - 1) * bSpline(n-1, u-1.0);
        }
    }
    
    /** Gets scaled fractional coordinate u as per Eq. 3.1 */
    public def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
        return new Vector3d(edgeReciprocals(0).mul(gridSize(0)).dot(r), edgeReciprocals(1).mul(gridSize(1)).dot(r), edgeReciprocals(2).mul(gridSize(2)).dot(r));
    }

    /**
     * Gets the volume V of the unit cell.
     */
    public def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }
}
