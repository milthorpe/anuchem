package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.fft.Distributed3dFft;
import au.edu.anu.util.Timer;
import x10.util.GrowableRail;
import x10.util.HashMap;
import x10.util.Pair;

/**
 * This class implements a Smooth Particle Mesh Ewald method to calculate
 * the potential of a system of charged particles.
 * TODO dipole correction for forces as described in Lambert, Darden & Board and Deem et al.
 * @see Essmann et al. "A Smooth Particle Mesh Ewald method", J. Comp. Phys. 101,
 * pp.8577-8593 (1995) DOI: 10.1063/1.470117
 */
public class PME {
    // TODO enum - XTENLANG-1118
    public const TIMER_INDEX_TOTAL : Int = 0;
    public const TIMER_INDEX_DIVIDE : Int = 1;
    public const TIMER_INDEX_DIRECT : Int = 2;
    public const TIMER_INDEX_SELF : Int = 3;
    public const TIMER_INDEX_GRIDCHARGES : Int = 4;
    public const TIMER_INDEX_INVFFT : Int = 5;
    public const TIMER_INDEX_THETARECCONVQ : Int = 6;
    public const TIMER_INDEX_RECIPROCAL : Int = 7;
    public const TIMER_INDEX_SETUP : Int = 8;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(9);

    /** The number of grid lines in each dimension of the simulation unit cell. */
    private global val gridSize : ValRail[Int](3);

    /** Double representations of the various grid dimensions */
    private global val K1 : Double;
    private global val K2 : Double; 
    private global val K3 : Double;

    /** The edges of the unit cell. */
    private global val edges : ValRail[Vector3d](3);
    private global val edgeLengths : ValRail[Double](3);

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

    /** 
     * Translation vectors for neighbouring unit cells 
     * (the 26 cells surrounding the origin cell)
     * Dimensions are:
     * 0: a "virtual" dimension used to replicate the data across all places
     *   // TODO should be done as a global array XTENLANG-787
     * 1: x translation (difference between x-coordinate of sub-cells
     * 2: y translation
     * 3: z translation
     */
    private global val imageTranslations : DistArray[Vector3d](4){rect};

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : DistArray[ValRail[MMAtom]](1);

    private global val B : DistArray[Double]{self.dist==gridDist};
    private global val C : DistArray[Double]{self.dist==gridDist};
    private global val BdotC : DistArray[Double]{self.dist==gridDist};

    private global val Q : DistArray[Complex]{self.dist==gridDist};

    /** 
     * An array of box divisions within the unit cell, with a side length
     * equal to the direct sum cutoff distance.  (N.B. if the unit cell side
     * length is not an exact multiple of the cutoff distance, the last box
     * in each dimension will be smaller than the cutoff distance, resulting
     * in anisotropy in the direct potential.)
     * Direct sums are only calculated between particles in the same box and
     * the 26 neighbouring boxes.
     * Dimensions of the array region are (x,y,z)
     * TODO assumes cubic unit cell
     */
    private global val subCells : PeriodicDistArray[ValRail[MMAtom]](3);
    /** The number of sub-cells per side of the unit cell. */
    private global val numSubCells : Int;

    /** 
     * A cache of "packed" representations of atoms from subcells
     * stored at other places.  This is used to prefetch atom data
     * for direct energy calculation.
     */
    private global val packedAtomsCache : DistArray[HashMap[Point(3), ValRail[MMAtom.PackedRepresentation]]](1);

    // TODO should be shared local to calculateEnergy() - XTENLANG-404
    private var directEnergy : Double = 0.0;
    private var selfEnergy: Double = 0.0;
    private var correctionEnergy : Double = 0.0; // TODO masklist
    private var reciprocalEnergy : Double = 0.0;

    /**
     * Creates a new particle mesh Ewald method.
     * @param edges the edge vectors of the unit cell
     * @param gridSize the number of grid lines in each dimension of the unit cell
     * @param atoms the atoms in the unit cell
     * @param splineOrder the order n of B-spline interpolationb
     * @param beta the Ewald coefficient beta
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(edges : ValRail[Vector3d](3),
            gridSize : ValRail[Int](3),
            atoms: DistArray[ValRail[MMAtom]](1),
            splineOrder : Int,
            beta : Double,
            cutoff : Double) {
        this.gridSize = gridSize;
        K1 = gridSize(0) as Double;
        K2 = gridSize(1) as Double;
        K3 = gridSize(2) as Double;
        this.edges = edges;
        this.edgeLengths = ValRail.make[Double](3, (i : Int) => edges(i).length());
        this.edgeReciprocals = ValRail.make[Vector3d](3, (i : Int) => edges(i).inverse());

        this.atoms = atoms;
        val r = Region.makeRectangular(0, gridSize(0)-1);
        gridRegion = (r * [0..(gridSize(1)-1)] * [0..(gridSize(2)-1)]) as Region(3);
        gridDist = Dist.makeBlockBlock(gridRegion, 0, 1);
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        val imageTranslationRegion = [0..Place.MAX_PLACES-1,-1..1,-1..1,-1..1] as Region(4){rect};
        this.imageTranslations = DistArray.make[Vector3d](Dist.makeBlock(imageTranslationRegion, 0), ((place,i,j,k) : Point(4)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k)));

        if (edgeLengths(0) % cutoff != 0.0) {
            Console.ERR.println("warning: edge length " + edgeLengths(0) + " is not an exact multiple of cutoff " + cutoff);
        }
        val numSubCells = Math.ceil(edgeLengths(0) / cutoff) as Int;
        val subCellRegion = [0..numSubCells-1,0..numSubCells-1,0..numSubCells-1] as Region(3){rect};
        val subCells = PeriodicDistArray.make[ValRail[MMAtom]](Dist.makeBlockBlock(subCellRegion, 0, 1));
        Console.OUT.println("subCells dist = " + subCells.dist);
        this.subCells = subCells;
        this.numSubCells = numSubCells;

        packedAtomsCache = DistArray.make[HashMap[Point(3), ValRail[MMAtom.PackedRepresentation]]](Dist.makeUnique(), (Point)=> new HashMap[Point(3), ValRail[MMAtom.PackedRepresentation]]());

        Console.OUT.println("gridDist = " + gridDist);

        Q = DistArray.make[Complex](gridDist);

        B = DistArray.make[Double](gridDist);
        C = DistArray.make[Double](gridDist);
        BdotC = DistArray.make[Double](gridDist);
    }

    /**
     * This method sets up the B, C and BdotC arrays, which only need be done once
     * in the simulation (as it is dependent on the grid size).
     * TODO this should be done in the constructor, but can't be done without big 
     * memory leaks due to proto rules.
     */
    public def setup() {
        timer.start(TIMER_INDEX_SETUP);
        initBArray();
        initCArray();
        finish ateach (p in BdotC) {
            BdotC(p) = B(p) * C(p);
        }
        timer.start(TIMER_INDEX_SETUP);
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        divideAtomsIntoSubCells();
        directEnergy = getDirectEnergy();
        selfEnergy = getSelfEnergy();
        gridCharges();

        timer.start(TIMER_INDEX_INVFFT);
        val temp = DistArray.make[Complex](gridDist); 
        val Qinv = DistArray.make[Complex](gridDist);
        new Distributed3dFft(gridSize(0), Q, Qinv, temp).doFFT3d(false);
        timer.stop(TIMER_INDEX_INVFFT);

        timer.start(TIMER_INDEX_THETARECCONVQ);
        // create F^-1(thetaRecConvQ)
        val thetaRecConvQ = DistArray.make[Complex](gridDist, (p : Point) => BdotC(p) * Qinv(p));

        // and do inverse FFT
        new Distributed3dFft(gridSize(0), thetaRecConvQ, thetaRecConvQ, temp).doFFT3d(true);
        timer.stop(TIMER_INDEX_THETARECCONVQ);
        reciprocalEnergy = getReciprocalEnergy(thetaRecConvQ);

        Console.OUT.println("directEnergy = " + directEnergy);
        Console.OUT.println("selfEnergy = " + selfEnergy);
        //Console.OUT.println("correctionEnergy = " + correctionEnergy);
        Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val totalEnergy = directEnergy + reciprocalEnergy + (correctionEnergy + selfEnergy);

        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    private def divideAtomsIntoSubCells() {
        timer.start(TIMER_INDEX_DIVIDE);
        val subCellsTemp = DistArray.make[GrowableRail[MMAtom]](subCells.dist, (Point) => new GrowableRail[MMAtom]());
        finish ateach (p1 in atoms) {
            val localAtoms = atoms(p1);
            foreach ((i) in 0..localAtoms.length-1) {
                val atom = localAtoms(i);
                val charge = atom.charge;
                val centre = atom.centre;
                val index = getSubCellIndex(atom);
                at (subCellsTemp.dist(index)) {
                    val remoteAtom = new MMAtom(centre, charge);
                    atomic{subCellsTemp(index).add(remoteAtom);}
                }
            }
        }
        finish ateach (p in subCells.dist) {
            subCells(p) = subCellsTemp(p).toValRail();
        }
        timer.stop(TIMER_INDEX_DIVIDE);
    }

    /**
     * Gets the index of the sub-cell for this atom within
     * the grid of sub-cells for direct sum calculation.
     */
    private global safe def getSubCellIndex(atom : MMAtom) : Point(3) {
        // TODO assumes cubic unit cell
        val r = Vector3d(atom.centre);
        val i = (r.i / cutoff) as Int;
        val j = (r.j / cutoff) as Int;
        val k = (r.k / cutoff) as Int;
        return Point.make(i, j, k);
    }

    public def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        finish ateach (p in subCells.dist) {
            var myDirectEnergy : Double = 0.0;
            val thisCell = subCells(p);
            for (var i : Int = p(0)-1; i<=p(0); i++) {
                var n1 : Int = 0;
                if (i < 0) {
                    n1 = -1;
                } // can't have (i > numSubCells+1)
                for (var j : Int = p(1)-1; j<=p(1)+1; j++) {
                    var n2 : Int = 0;
                    if (j < 0) {
                        n2 = -1;
                    } else if (j > numSubCells-1) {
                        n2 = 1;
                    }
                    for (var k : Int = p(2)-1; k<=p(2)+1; k++) {
                        var n3 : Int = 0;
                        if (k < 0) {
                            n3 = -1;
                        } else if (k > numSubCells-1) {
                            n3 = 1;
                        }
                        // interact with "left half" of other boxes i.e. only boxes with i<=p(0)
                        if (i < p(0) || (i == p(0) && j < p(1)) || (i == p(0) && j == p(1) && k < p(2))) {
                            val translation = imageTranslations(here.id,n1,n2,n3);
                            val otherSubCellLocation = subCells.periodicDist(i,j,k);
                            if (otherSubCellLocation == here) {
                                val otherCell = subCells(i,j,k);
                                for ((otherAtom) in 0..otherCell.length()-1) {
                                    val imageLoc = otherCell(otherAtom).centre + translation;
                                    for ((thisAtom) in 0..thisCell.length()-1) {
                                        val r = thisCell(thisAtom).centre.distance(imageLoc);
                                        if (r < cutoff) {
                                            val chargeProduct = thisCell(thisAtom).charge * otherCell(otherAtom).charge;
                                            val imageDirectComponent = chargeProduct * Math.erfc(beta * r) / r;
                                            myDirectEnergy += imageDirectComponent;
                                        }
                                    }
                                }
                            } else {
                                // other subcell is remote; need to transfer packed atoms
                                var otherCellPacked : ValRail[MMAtom.PackedRepresentation] = null;
                                // cache remote atoms
                                val otherCellIndex = Point.make(i,j,k);
                                atomic {otherCellPacked = packedAtomsCache(here.id).getOrElse(otherCellIndex, null);}
                                if (otherCellPacked == null) {
                                    otherCellPacked = at (subCells.periodicDist(otherCellIndex)) {getPackedAtomsForSubCell(otherCellIndex)};
                                    atomic {packedAtomsCache(here.id).put(otherCellIndex, otherCellPacked);}
                                }
                                for ((otherAtom) in 0..otherCellPacked.length()-1) {
                                    val imageLoc = otherCellPacked(otherAtom).centre + translation;
                                    for ((thisAtom) in 0..thisCell.length()-1) {
                                        val r = thisCell(thisAtom).centre.distance(imageLoc);
                                        if (r < cutoff) {
                                            val chargeProduct = thisCell(thisAtom).charge * otherCellPacked(otherAtom).charge;
                                            val imageDirectComponent = chargeProduct * Math.erfc(beta * r) / r;
                                            myDirectEnergy += imageDirectComponent;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // atoms in same cell
            for ((i) in 0..thisCell.length()-1) {
                for ((j) in 0..i-1) {
                    val rjri = thisCell(j).centre - thisCell(i).centre;
                    val r = rjri.length();
                    if (r < cutoff) {
                        val directComponent = thisCell(i).charge * thisCell(j).charge * Math.erfc(beta * r) / r;
                        myDirectEnergy += directComponent;
                    }
                }
            }

            // TODO this is slow because of lack of optimized atomic - XTENLANG-321
            val myDirectEnergyFinal = myDirectEnergy;
            at(this) {
               atomic directEnergy += myDirectEnergyFinal;
            }
        }
        
        directEnergy = directEnergy;
        timer.stop(TIMER_INDEX_DIRECT);
        return directEnergy;
    }

    /*
     * Returns atom charges and coordinates for a sub-cell, in packed representation
     */
    private global safe def getPackedAtomsForSubCell(subCellIndex : Point(3)) : ValRail[MMAtom.PackedRepresentation] {
        val subCell = subCells(subCellIndex);
        return ValRail.make[MMAtom.PackedRepresentation](subCell.length(), (i : Int) => subCell(i).getPackedRepresentation());
    }

    /**
     * Returns the self energy as defined in Eq. 2.5, which is the
     * contribution of the interaction in reciprocal space of each
     * atom with itself.  This is subtracted from the final energy.
     */
    public def getSelfEnergy() : Double {
        timer.start(TIMER_INDEX_SELF);
        finish ateach (p in subCells.dist) {
            val thisCell = subCells(p);
            var mySelfEnergy : Double = 0.0;
            for ((thisAtom) in 0..thisCell.length()-1) {
                mySelfEnergy += thisCell(thisAtom).charge * thisCell(thisAtom).charge;
            }
            val mySelfEnergyFinal = mySelfEnergy;
            at(this) {
                atomic selfEnergy += mySelfEnergyFinal;
            }
        }
        selfEnergy = -beta / Math.sqrt(Math.PI) * selfEnergy;
        timer.stop(TIMER_INDEX_SELF);
        return selfEnergy;
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     */
    public def gridCharges() {
        timer.start(TIMER_INDEX_GRIDCHARGES);
        finish ateach ((place1) in Dist.makeUnique(gridDist.places())) {
            val myQ = new PeriodicArray[Double](gridRegion);
            for (p in subCells | here) {
                val thisCell = subCells(p);
                for (atom in thisCell) {
                    val q = atom.charge;
                    val u = getScaledFractionalCoordinates(atom.centre); // TODO general non-cubic
                    val u1c = Math.ceil(u.i) as Int;
                    val u2c = Math.ceil(u.j) as Int;
                    val u3c = Math.ceil(u.k) as Int;
                    for ((i) in 1..splineOrder) {
                        val k1 = (u1c - i);
                        for ((j) in 1..splineOrder) {
                            val k2 = (u2c - j);
                            for ((k) in 1..splineOrder) {
                                val k3 = (u3c - k);
                                val gridPointContribution = q
                                             * bSpline(splineOrder, u.i - k1)
                                             * bSpline(splineOrder, u.j - k2)
                                             * bSpline(splineOrder, u.k - k3);
                                myQ(k1,k2,k3) = myQ(k1,k2,k3) + gridPointContribution;
                            }
                        }
                    }
                }
            }
            // scatter myQ and accumulate to distributed Q
            finish foreach (place2 in gridDist.places()) {
                val place1ContributionDist = Q.dist | place2;
                val place1Contribution = Rail.make[Double](place1ContributionDist.region.size());
                var i : Int = 0;
                for (p in place1ContributionDist) {
                    place1Contribution(i++) = myQ(p);
                }
                val place1ContributionToTransfer = place1Contribution as ValRail[Double];
                at (place2) {
                    atomic {
                        var j : Int = 0;
                        for (p in Q | here) {
                            Q(p) = Q(p) + place1ContributionToTransfer(j++);
                        }
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_GRIDCHARGES);
    }

    /**
     * @param Q the gridded charge array as defined in Eq. 4.6
     * @param thetaRecConvQ the reciprocal pair potential as defined in eq. 4.7
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(thetaRecConvQ : DistArray[Complex]{self.dist==gridDist}) {
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
     * Initialises the array B as defined by Eq 4.8 and 4.4
     * TODO should be able to construct array as local and return, but can't due to GC limitations
     */
    public def initBArray() {
        finish ateach ((m1,m2,m3) in B) {
            val m1D = m1 as Double;
            var sumK1 : Complex = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK1 = sumK1 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m1D * k / K1 * Complex.I);
            }
            val b1 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m1D / K1 * Complex.I) / sumK1).abs();

            val m2D = m2 as Double;
            var sumK2 : Complex = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK2 = sumK2 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m2D * k / K2 * Complex.I);
            }
            val b2 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m2D / K2 * Complex.I) / sumK2).abs();
                
            val m3D = m3 as Double;
            var sumK3 : Complex = Complex.ZERO;
            for ((k) in 0..(splineOrder-2)) {
                sumK3 = sumK3 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
            }
            val b3 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m3D / K3 * Complex.I) / sumK3).abs();

            B(m1,m2,m3) = b1 * b1 * b2 * b2 * b3 * b3;
            //Console.OUT.println("B(" + m1 + "," + m2 + "," + m3 + ") = " + B(m1,m2,m3));
        }
    }

    /**
     * Initialises the array C as defined by Eq 3.9
     * TODO should be able to construct array as local and return, but can't due to GC limitations
     */
    public def initCArray() {
        val V = getVolume();
        //Console.OUT.println("V = " + V);
        finish ateach ((m1,m2,m3) in C) {
            val m1prime = m1 <= K1/2 ? m1 : m1 - K1;
            val m2prime = m2 <= K2/2 ? m2 : m2 - K2;
            val m3prime = m3 <= K3/2 ? m3 : m3 - K3;
            val mVec = edgeReciprocals(0).mul(m1prime).add(edgeReciprocals(1).mul(m2prime)).add(edgeReciprocals(2).mul(m3prime));
            val mSquared = mVec.dot(mVec);
            C(m1,m2,m3) = Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)) / (mSquared * Math.PI * V);
        }
        at (C.dist(0,0,0)) {
            C(0,0,0) = 0.0;
        }
    }

    /* 
     * Gets the nth order B-spline M_n(u) as per Eq. 4.1
     */
    public final static safe def bSpline(n : Int, u : Double) : Double {
        if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 4) {
            return bSpline4(u);
        } else if (n == 2) {
            return 1.0 - Math.abs(u - 1.0);
        } else {
            return u / (n - 1) * bSpline(n-1, u) + (n - u) / (n - 1) * bSpline(n-1, u-1.0);
        }
    }

    /* 
     * Gets the 4th order B-spline M_4(u) as per Eq. 4.1
     */
    private final static safe def bSpline4(u : Double) : Double {
        if (u <= 0.0 || u >= 4) {
            return 0.0;
        } else {
            return u / 3 * bSpline3(u) + (4 - u) / 3 * bSpline3(u-1.0);
        }
    }

    /* 
     * Gets the 3rd order B-spline M_3(u) as per Eq. 4.1
     */
    private final static safe def bSpline3(u : Double) : Double {
        if (u <= 0.0 || u >= 3) {
            return 0.0;
        } else {
            return u / 2 * bSpline2(u) + (3 - u) / 2 * bSpline2(u-1.0);
        }
    }

    /* 
     * Gets the 2nd order B-spline M_2(u) as per Eq. 4.1
     */
    private final static safe def bSpline2(u : Double) : Double {
        if (u <= 0.0 || u >= 2) {
            return 0.0;
        } else {
            return 1.0 - Math.abs(u - 1.0);
        }
    }

    /** Gets scaled fractional coordinate u as per Eq. 3.1 - cubic only */
    public global safe def getScaledFractionalCoordinates(r : Point3d) : Vector3d {
        return Vector3d(edgeReciprocals(0).i * K1 * r.i, edgeReciprocals(1).j * K2 * r.j, edgeReciprocals(2).k * K3 * r.k);
    }
    
    /** Gets scaled fractional coordinate u as per Eq. 3.1 - general rectangular */
    public global safe def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
        // this method allows general non-rectangular cells
        return Vector3d(edgeReciprocals(0).mul(K1).dot(r), edgeReciprocals(1).mul(K2).dot(r), edgeReciprocals(2).mul(K3).dot(r));
    }

    /**
     * Gets the volume V of the unit cell.
     */
    private global safe def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }
}
