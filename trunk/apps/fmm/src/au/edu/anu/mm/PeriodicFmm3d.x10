/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.util.*;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * This subclass of Fmm3d extends the base FMM with periodic boundary
 * conditions.  In addition to the interactions within the unit cell, 
 * the  unit cell interacts with 3^k * 3^k * 3^k copies of itself in 
 * concentric shells of increasingly coarse-grained aggregate cells.
 *
 * @see Lambert, Darden & Board (1996). "A Multipole-Based Algorithm 
 *  for Efficient Calculation of Forces and Potentials in Macroscopic 
 *  Periodic Assemblies of Particles". J Comp Phys 126 274-285
 *
 * @see Kudin & Scuseria (1998). "A fast multipole method for Periodic 
 * systems with arbitrary unit cell geometries". 
 * Chem. Phys. Letters 283 61-68
 * @author milthorpe
 */
public class PeriodicFmm3d extends Fmm3d {
    /** The number of concentric shells of copies of the unit cell. */
    public global val numShells : Int;

    /** All boxes in the octree division of space. 
     * ValRail has numLevels elements, for levels [1..numLevels]
     * where boxes(n) = boxes at level n+1
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    global val boxes : ValRail[PeriodicDistArray[FmmBox](3)];

    global val lowestLevelBoxes : PeriodicDistArray[FmmBox](3);

    /**
     * An array of locally essential trees (LETs), one for each place.
     */
    protected global val locallyEssentialTrees : DistArray[LocallyEssentialTree](1);

    public const TIMER_INDEX_MACROSCOPIC : Int = 7;
    /** 
     * A multi-timer for the several segments of a single getEnergy 
     * invocation, indexed by the constants above and in the superclass. 
     */
    public val timer = new Timer(8);

    /** A region representing a cube of 3x3x3 boxes, used for constructing macroscopic multipoles. */
    global val threeCube : Region(3) = [-1..1,-1..1,-1..1] as Region(3);

    /** A region representing a cube of 9x9x9 boxes, used for interacting with macroscopic multipoles. */
    global val nineCube : Region(3) = [-4..4,-4..4,-4..4] as Region(3);

    /** The net dipole moment of the unit cell. */
    var dipole : Vector3d = Vector3d.NULL;

    /**
     * Initialises a periodic fast multipole method electrostatics 
     * calculation for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     * @param numShells the number of concentric shells of copies of the unit cell
     */
    public def this(density : Double, 
                    numTerms : Int,
                    topLeftFront : Point3d,
                    size : Double,  
                    numAtoms : Int,
                    atoms: DistArray[ValRail[MMAtom]](1),
                    numShells : Int) {
        // Periodic FMM always uses ws = 1
        // TODO is it possible to formulate for well-spaced > 1?
        super(density, numTerms, 1, topLeftFront, size, numAtoms, atoms);
        timer.start(TIMER_INDEX_TREE);
        val boxes = constructTree();
        this.boxes = boxes;
        this.lowestLevelBoxes = boxes(numLevels);
        this.numShells = numShells;
        this.locallyEssentialTrees = createLocallyEssentialTrees();
        assignAtomsToBoxes(atoms, boxes(numLevels));
        timer.stop(TIMER_INDEX_TREE);
    }

    public global def getTopLevel() : Int = 0;
    
    public def calculateEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);
        multipoleLowestLevel();
        // direct energy is independent of all subsequent steps of FMM
        val directEnergy = future {getDirectEnergy()};
        combineMultipoles();
        combineMacroscopicExpansions();
        transformToLocal();
        farFieldEnergy = getFarFieldEnergy();
        val totalEnergy = directEnergy() + farFieldEnergy;
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    /** 
     * Generate and includes macroscopic expansions, for concentric rings
     * of aggregates of copies of the unit cell.
     */
    def combineMacroscopicExpansions() {
        timer.start(TIMER_INDEX_MACROSCOPIC);
        // TODO distributed impl.
        val macroMultipoles = DistArray.make[MultipoleExpansion](Dist.makeBlock([0..Place.MAX_PLACES-1,0..numShells],0));
        val macroLocalTranslations = DistArray.make[LocalExpansion](Dist.makeBlock([0..Place.MAX_PLACES-1,0..numShells],0));
        finish ateach ((p) : Point in Dist.makeUnique(Place.places)) {
            //Console.OUT.println("shells at " + p + " top level = " + getTopLevel());
            val topLevelBox = boxes(0)(0,0,0) as FmmBox!;
            // "shell 0" is the unit cell itself
            macroMultipoles(p, 0) = topLevelBox.multipoleExp;

            var macroTranslation : MultipoleExpansion! = new MultipoleExpansion(numTerms);

            // multipoles for shell 1
            for ((i,j,k) in threeCube) {
                val translationVector = Vector3d(i * size,
                                                 j * size,
                                                 k * size);
                val translation = MultipoleExpansion.getOlm(translationVector, numTerms);
                macroTranslation.add(translation);
            }
            macroMultipoles(p, 1) = new MultipoleExpansion(numTerms);
            macroMultipoles(p, 1).translateAndAddMultipole(macroTranslation, macroMultipoles(p, 0));
            //Console.OUT.println("final for 1 = " + macroMultipoles(p, 1));

            // locals for shell 1
            macroLocalTranslations(p,0) = new LocalExpansion(numTerms);
            for ((i,j,k) in nineCube) {
                if (Math.abs(i) > 1 || Math.abs(j) > 1 || Math.abs(k) > 1) {
                    // inner 27 boxes done at a lower level
                    val translationVector = Vector3d(i * size,
                                                     j * size,
                                                     k * size);
                    val transform = LocalExpansion.getMlm(translationVector, numTerms) as LocalExpansion!;
                    macroLocalTranslations(p,0).add(transform);
                }
            }
            macroLocalTranslations(p,1) = macroLocalTranslations(p,0).getMacroscopicParent();

            // remaining shells
            for (var shell: Int = 2; shell <= numShells; shell++) {
                macroTranslation = macroTranslation.getMacroscopicParent();
                macroMultipoles(p, shell) = new MultipoleExpansion(numTerms);
                macroMultipoles(p, shell).translateAndAddMultipole(macroTranslation, macroMultipoles(p, shell-1));
                //Console.OUT.println("final for " + shell + " = " + macroMultipoles(p, shell));
                macroLocalTranslations(p,shell) = macroLocalTranslations(p,shell-1).getMacroscopicParent();
            }

            // now transform and add macroscopic multipoles to local expansion for top level box
            for (var shell: Int = 0; shell <= numShells; shell++) {
                val localExpansion = macroLocalTranslations(p,shell) as LocalExpansion!;
                topLevelBox.localExp.transformAndAddToLocal(localExpansion, macroMultipoles(p, shell) as MultipoleExpansion!);
            }
            //Console.OUT.println("final for topLevel = " + topLevelBox.localExp);
        }
        timer.stop(TIMER_INDEX_MACROSCOPIC);
    }

    def assignAtomsToBoxes(atoms: DistArray[ValRail[MMAtom]](1), lowestLevelBoxes : PeriodicDistArray[FmmBox](3)) {
        //Console.OUT.println("assignAtomsToBoxes");
        finish ateach (p1 in atoms) {
            var myDipole : Vector3d = Vector3d.NULL;
            val localAtoms = atoms(p1);
            for ((i) in 0..localAtoms.length-1) {
                val atom = localAtoms(i) as MMAtom!;
                val charge = atom.charge;
                val offsetCentre = atom.centre + offset;
                myDipole = myDipole + Vector3d(offsetCentre) * charge;
                val boxIndex = getLowestLevelBoxIndex(offsetCentre);
                async(lowestLevelBoxes.dist(boxIndex)) {
                    val remoteAtom = new MMAtom(offsetCentre, charge);
                    val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
                    atomic {
                        leafBox.atoms.add(remoteAtom);
                    }
                }
            }
            val myDipoleFinal = myDipole;
            at (this) { atomic { dipole = dipole + myDipoleFinal; } }
        }

        // post-prune leaf boxes
        // TODO prune intermediate empty boxes as well
        finish ateach (boxIndex in lowestLevelBoxes) {
            val box = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
            if (box.atoms.length() == 0) {
                lowestLevelBoxes(boxIndex) = null;
            }
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        //Console.OUT.println("multipole lowest level");
        timer.start(TIMER_INDEX_MULTIPOLE);
        finish ateach (boxIndex in lowestLevelBoxes) {
            val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
            val boxLocation = leafBox.getCentre(size);
            for ((i) in 0..leafBox.atoms.length()-1) {
                val atom = leafBox.atoms(i);
                val atomLocation = leafBox.getCentre(size).vector(atom.centre);
                val atomExpansion = MultipoleExpansion.getOlm(atom.charge, atomLocation, numTerms);
                leafBox.multipoleExp.add(atomExpansion);
            }
        }

        cancelDipole(dipole);

        timer.stop(TIMER_INDEX_MULTIPOLE);
    }

    def addAtomToLowestLevelBoxAsync(boxIndex : Point(3), offsetCentre : Point3d, charge : Double) {
        async(lowestLevelBoxes.dist(boxIndex)) {
            val remoteAtom = new MMAtom(offsetCentre, charge);
            val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
            val boxLocation = leafBox.getCentre(size).vector(offsetCentre);
            val atomExpansion = MultipoleExpansion.getOlm(charge, boxLocation, numTerms);
            atomic {
                leafBox.atoms.add(remoteAtom);
                leafBox.multipoleExp.add(atomExpansion);
            }
        }
    }

    /** 
     * Add fictious charges to the corners of the unit cell 
     * to cancel the dipole moment.
     * @see Kudin & Scuseria, section 2.3
     */
    def cancelDipole(dipole : Vector3d) : Vector3d {
        //Console.OUT.println("dipole = " + dipole);
        var newDipole : Vector3d = dipole;
        finish {
            val p1 = Point3d(size, 0.0, 0.0) + offset;
            val q1 = - dipole.i / size;
            addAtomToLowestLevelBoxAsync(Point.make(lowestLevelDim-1, 0, 0), p1, q1);
            newDipole = newDipole + Vector3d(p1) * q1;

            val p2 = Point3d(0.0, size, 0.0) + offset;
            val q2 = - dipole.j / size;
            addAtomToLowestLevelBoxAsync(Point.make(0, lowestLevelDim-1, 0), p2, q2);
            newDipole = newDipole + Vector3d(p2) * q2;


            val p3 = Point3d(0.0, 0.0, size) + offset;
            val q3 = - dipole.k / size;
            addAtomToLowestLevelBoxAsync(Point.make(0, 0, lowestLevelDim-1), p3, q3);
            newDipole = newDipole + Vector3d(p3) * q3;


            val p0 = Point3d(0.0, 0.0, 0.0)  + offset;
            val q0 = -(q1 + q2 + q3);
            addAtomToLowestLevelBoxAsync(Point.make(0, 0, 0),                p0, q0);
            newDipole = newDipole + Vector3d(p0) * q0;
/*
            Console.OUT.println(q1 + " at " + p1);
            Console.OUT.println(q2 + " at " + p2);
            Console.OUT.println(q3 + " at " + p3);
            Console.OUT.println(q0 + " at " + p0);
*/
        }

        //Console.OUT.println("after cancelling, dipole = " + newDipole);
        return newDipole; 
    }

    /**
     * Starting at the top level, for each box, transform multipole
     * expansions for all well-separated boxes (for which the parent
     * box is not also well-separated) to local expansions for this 
     * box.  Translate the local expansion for each parent down to all
     * non-empty child boxes.
     * Note: this periodic version includes boxes in 8 periodic images
     * of the unit cell.
     */
    def transformToLocal() {
        timer.start(TIMER_INDEX_TRANSFORM);

        // start the prefetch of all multipoles required at each place
        prefetchMultipoles();
        
        for ((thisLevel) in (getTopLevel()+1)..numLevels) {
            val dim = size / Math.pow2(thisLevel);
            //Console.OUT.println("transform level " + thisLevel);
            val thisLevelBoxes = boxes(thisLevel);
            finish ateach (p1 in Dist.makeUnique(thisLevelBoxes.dist.places())) {
                val thisLevelMultipoleCopies = (locallyEssentialTrees(here.id) as LocallyEssentialTree!).multipoleCopies(thisLevel);
                foreach ((x1,y1,z1) in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(x1,y1,z1) as FmmBox!;
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for ((x2,y2,z2) in vList) {
                            // force on the multipole value
                            val box2MultipoleExp = thisLevelMultipoleCopies(x2,y2,z2)() as MultipoleExpansion!;
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(thisLevel,x2,y2,z2);
                                val translateP = Point.make([here.id, translation(0), translation(1), translation(2), translation(3)]);
                                val transform21 = multipoleTransforms(translateP) as LocalExpansion!;
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                        val box1ParentExp = box1.parent.getLocalExpansionLocalCopy(numTerms);
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, box1.x%2, box1.y%2, box1.z%2])) as MultipoleExpansion!;
                        box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_TRANSFORM);
    }

    private def constructTree() : ValRail[PeriodicDistArray[FmmBox](3)] {
        val boxesTemp = Rail.make[PeriodicDistArray[FmmBox](3)](numLevels+1);
        for ((thisLevel) in getTopLevel()..numLevels) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(3) = [0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = Dist.makeBlock(thisLevelRegion, 0);
            boxesTemp(thisLevel) = PeriodicDistArray.make[FmmBox](thisLevelDist);
            //Console.OUT.println("level " + thisLevel + " dist: " + thisLevelDist);
        }
        val boxesValRail = boxesTemp as ValRail[PeriodicDistArray[FmmBox](3)];

        for ((thisLevel) in getTopLevel()..numLevels-1) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelBoxes = boxesValRail(thisLevel);
            finish ateach ((x,y,z) in thisLevelBoxes) {
                val box = new FmmBox(thisLevel, x, y, z, numTerms, getParentForChild(boxesValRail, thisLevel, x,y,z));
                createVList(box);
                thisLevelBoxes(x,y,z) = box;
            }
        }

        val lowestLevelBoxes = boxesValRail(numLevels);
        finish ateach ((x,y,z) in lowestLevelBoxes) {
            val box = new FmmLeafBox(numLevels, x, y, z, numTerms, getParentForChild(boxesValRail, numLevels, x,y,z));
            createUList(box);
            createVList(box);
            lowestLevelBoxes(x,y,z) = box;
        }

        return boxesValRail;
    }

    /**
     * Creates the U-list of <code>box</code>.
     * The U-list consists of all leaf boxes not well-separated from <code>box</code>.
     */
    private global def createUList(box : FmmLeafBox!) {
        // interact with "left half" of uList i.e. only boxes with x<=box.x
        val uList = new GrowableRail[Point(3)]();
        for ((x) in box.x-ws..box.x) {
            for ((y) in box.y-ws..box.y+ws) {
                for ((z) in box.z-ws..box.z+ws) {
                    if (x < box.x || (x == box.x && y < box.y) || (x == box.x && y == box.y && z < box.z)) {
                        uList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setUList(uList.toValRail());
    }

    /**
     * Creates the V-list of <code>box</code>.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent of <code>box</code>.
     */
    private global def createVList(box : FmmBox!) {
        val xOffset = box.x%2 == 1 ? -1 : 0;
        val yOffset = box.y%2 == 1 ? -1 : 0;
        val zOffset = box.z%2 == 1 ? -1 : 0;
        val vList = new GrowableRail[Point(3)]();
        for ((x) in box.x-2*ws+xOffset..box.x+2*ws+1+xOffset) {
            for ((y) in box.y-2*ws+yOffset..box.y+2*ws+1+yOffset) {
                for ((z) in box.z-2*ws+zOffset..box.z+2*ws+1+zOffset) {
                    if (box.wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setVList(vList.toValRail());
    }

    /**
     * "Pre-fetches" the multipole copies for the V-List of the locally
     * essential tree at each place, using futures.
     */
    protected def prefetchMultipoles() {
        finish ateach (p1 in locallyEssentialTrees) {
            val myLET = locallyEssentialTrees(p1) as LocallyEssentialTree!;
            val combinedVList = myLET.combinedVList;
            for ((level) in getTopLevel()..numLevels) {
                val thisLevelCopies = myLET.multipoleCopies(level);
                val thisLevelBoxes = boxes(level);
                for ((x,y,z) in combinedVList(level)) {
                    if (thisLevelCopies(x,y,z) == null) {
                        thisLevelCopies(x,y,z) = future {getMultipoleExpansionLocalCopy(thisLevelBoxes,x,y,z)};
                    }
                }
            }
        }
    }

    /**
     * Gets sum of direct (pairwise) energy for all pairs of atoms
     * in non-well-separated boxes. This operations requires only
     * that atoms have already been assigned to boxes, and so can 
     * be done in parallel with other steps of the algorithm.
     */
    def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);

        finish ateach (p1 in locallyEssentialTrees) {
            // start the prefetch of all atoms required for direct calculations at each place
            val myLET = locallyEssentialTrees(p1) as LocallyEssentialTree!;
            val myCombinedUList = myLET.combinedUList;
            val packedAtoms = myLET.packedAtoms;
            for ((x,y,z) in myCombinedUList) {
                myLET.packedAtoms(x,y,z) = future {at (lowestLevelBoxes.periodicDist(x,y,z)) {getPackedAtomsForBox(x, y, z)}};
            }

            foreach ((x1,y1,z1) in lowestLevelBoxes | here) {
                val box1 = lowestLevelBoxes(x1,y1,z1) as FmmLeafBox!;
                if (box1 != null) {
                    // TODO use shared var - XTENLANG-404
                    var thisBoxEnergy : Double = 0.0;
                    val length = box1.atoms.length();
                    for ((atomIndex1) in 0..length-1) {
                        val atom1 = box1.atoms(atomIndex1);

                        // direct calculation with all atoms in same box
                        for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                            val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                            val pairEnergy = atom1.charge * sameBoxAtom.charge / atom1.centre.distance(sameBoxAtom.centre);
                            thisBoxEnergy += pairEnergy;
                        }
                    }

                    // direct calculation with all atoms in non-well-separated boxes
                    val uList = box1.getUList();
                    for ((x2,y2,z2) in uList) {
                        // here we force on the packed atoms for which a future was previously issued
                        val boxAtoms = packedAtoms(x2,y2,z2)();
                        if (boxAtoms != null) {
                            // is box from a neighbouring image?
                            val translation = getTranslation(x2, y2, z2);
                            for (var i : Int = 0; i < boxAtoms.length(); i++) {
                                val atom2Packed = boxAtoms(i);
                                val imageCentre = atom2Packed.centre + translation;
                                for ((atomIndex1) in 0..length-1) {
                                    val atom1 = box1.atoms(atomIndex1);
                                    if (atom1.centre != imageCentre) {
                                    // necessary to check to avoid interactions between fictious charges
                                        thisBoxEnergy += atom1.charge * atom2Packed.charge / atom1.centre.distance(imageCentre);
                                    }
                                }
                            }
                        }
                    }
                    val thisBoxEnergyFinal = thisBoxEnergy;
                    async (this) {atomic {directEnergy += thisBoxEnergyFinal;}}
                }
            }
        }
        timer.stop(TIMER_INDEX_DIRECT);

        return directEnergy;
    }

    /**
     * Gets the atom centre translation vector due to a lowest-level box 
     * being in a neighbouring image, rather than the central unit cell.
     */
    def getTranslation(x:  Int, y : Int, z : Int) : Vector3d {
        var translationX : Double = 0.0;
        if (x >= lowestLevelDim) {
            translationX = size;
        } else if (x < 0) {
            translationX = -size;
        }

        var translationY : Double = 0;
        if (y >= lowestLevelDim) {
            translationY = size;
        } else if (y < 0) {
            translationY = -size;
        }

        var translationZ : Double = 0;
        if (z >= lowestLevelDim) {
            translationZ = size;
        } else if (z < 0) {
            translationZ = -size;
        }
        return Vector3d(translationX, translationY, translationZ);
    }

    //
    // TODO all after this point shouldn't be necessary if DistArray, PeriodicDistArray can inherit from same class
    //

    private global def getParentForChild(boxes : ValRail[PeriodicDistArray[FmmBox](3)], level : Int, x : Int, y : Int, z : Int) : FmmBox {
        if (level == getTopLevel())
            return null;
        return (at (boxes(level-1).dist(x/2, y/2, z/2)) {boxes(level-1)(x/2, y/2, z/2)});
    }

    /**
     * Creates the locally essential tree at each place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private def createLocallyEssentialTrees() : DistArray[LocallyEssentialTree](1) {
        val locallyEssentialTrees = DistArray.make[LocallyEssentialTree](Dist.makeUnique(), (Point)=> null);
        finish ateach ((p1) in locallyEssentialTrees) {
            val uMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
            val uMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
            val combinedUSet = new HashSet[Point(3)]();
            for ((x,y,z) in lowestLevelBoxes | here) {
                val box1 = lowestLevelBoxes(x,y,z) as FmmLeafBox!;
                if (box1 != null) {
                    val uList = box1.getUList();
                    for (boxIndex2 in uList) {
                        if (combinedUSet.add(boxIndex2)) {
                            for ((i) in 0..2) {
                                uMin(i) = Math.min(uMin(i), boxIndex2(i));
                                uMax(i) = Math.max(uMax(i), boxIndex2(i));        
                            }
                        }
                    }
                }
            }
            val uListMin = uMin as ValRail[Int](3);
            val uListMax = uMax as ValRail[Int](3);
            val combinedUList = combinedUSet.toValRail();

            val combinedVList = Rail.make[ValRail[Point(3)]](numLevels+1);
            val vListMin = Rail.make[ValRail[Int](3)](numLevels+1);
            val vListMax = Rail.make[ValRail[Int](3)](numLevels+1);
            for ((thisLevel) in getTopLevel()..numLevels) {
                val vMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
                val vMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
                //Console.OUT.println("create combined V-list for level " + thisLevel + " at " + here);
                val combinedVSet = new HashSet[Point(3)]();
                val thisLevelBoxes = boxes(thisLevel);
                for ((x,y,z) in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(x,y,z) as FmmBox!;
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for (boxIndex2 in vList) {
                            if (combinedVSet.add(boxIndex2)) {
                                for ((i) in 0..2) {
                                    vMin(i) = Math.min(vMin(i), boxIndex2(i));
                                    vMax(i) = Math.max(vMax(i), boxIndex2(i));        
                                }
                            }
                        }
                    }
                }
                //Console.OUT.println("done " + combinedVSet.size());
                vListMin(thisLevel) = vMin as ValRail[Int](3);
                vListMax(thisLevel) = vMax as ValRail[Int](3);
                combinedVList(thisLevel) = combinedVSet.toValRail();
            }
            locallyEssentialTrees(p1) = new LocallyEssentialTree(combinedUList,
                                                                 combinedVList as ValRail[ValRail[Point(3)]],
                                                                 uListMin,
                                                                 uListMax,
                                                                 vListMin as ValRail[ValRail[Int]],
                                                                 vListMax as ValRail[ValRail[Int]]);
        }
        return locallyEssentialTrees;
    }

    /**
     * TODO this is a workaround due to lack of Array copy facility - XTENLANG-787
     * @return a local copy at the current place of this box's multipole expansion
     */
    private global def getMultipoleExpansionLocalCopy(thisLevelBoxes : PeriodicDistArray[FmmBox](3), x : Int, y : Int, z : Int) : MultipoleExpansion! {
        val data = at (thisLevelBoxes.periodicDist(x,y,z)) {thisLevelBoxes(x,y,z) != null? Expansion.getData(numTerms, (thisLevelBoxes(x,y,z) as FmmBox!).multipoleExp) : null};
        if (data != null) {
            return new MultipoleExpansion(numTerms, data);
        } else {
            return null;
        }
    }

    private global def getPackedAtomsForBox(x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox!;
        if (box != null) {
            return box.getPackedAtoms();
        } else {
            return null;
        }
    }

    /** 
     * Starting at the bottom level, combines multipole expansions for <= 8 child
     * boxes into a single multipole expansion for the parent box.
     */
    def combineMultipoles() {
        timer.start(TIMER_INDEX_COMBINE);
        for (var level: Int = numLevels; level > getTopLevel(); level--) {
            val thisLevel = level;
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelBoxes = boxes(thisLevel);
            finish ateach (boxIndex in thisLevelBoxes) {
                val child = thisLevelBoxes(boxIndex) as FmmBox!;
                if (child != null) {
                    val childExp = child.multipoleExp;
                    val parent = child.parent;
                    at (parent) {
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, (child.x+1)%2, (child.y+1)%2, (child.z+1)%2])) as MultipoleExpansion!;
                        parent.multipoleExp.translateAndAddMultipole(shift, childExp);
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_COMBINE);
    }

    /**
     * Calculates the sum of far-field interactions for all atoms:
     * for each atom, calculates the interaction using the local expansion
     * for the enclosing box, which is the combined contribution of
     * all atoms in all well-separated boxes.
     */ 
    def getFarFieldEnergy() : Double {
        //Console.OUT.println("getFarFieldEnergy");
        timer.start(TIMER_INDEX_FARFIELD);
        finish ateach (boxIndex1 in lowestLevelBoxes) {
            val box1 = lowestLevelBoxes(boxIndex1) as FmmLeafBox!;
            if (box1 != null) {
                // TODO use shared var - XTENLANG-404
                var thisBoxEnergy : Double = 0.0;
                val length = box1.atoms.length();
                for ((atomIndex1) in 0..length-1) {
                    val atom1 = box1.atoms(atomIndex1);
                    val box1Centre = atom1.centre.vector(box1.getCentre(size));
                    val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                    thisBoxEnergy += farFieldEnergy;
                }
                val thisBoxEnergyFinal = thisBoxEnergy;
                async (this) {atomic {farFieldEnergy += thisBoxEnergyFinal;}}
            }
        }
        timer.stop(TIMER_INDEX_FARFIELD);

        return farFieldEnergy / 2.0;
    }

}

