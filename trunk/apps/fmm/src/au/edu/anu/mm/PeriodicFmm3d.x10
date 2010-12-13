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
    public val numShells : Int;

    /** All boxes in the octree division of space. 
     * Array has numLevels elements, for levels [1..numLevels]
     * where boxes(n) = boxes at level n+1
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    val boxes : Array[PeriodicDistArray[FmmBox]{rank==3,rect}](1){rail};

    val lowestLevelBoxes : PeriodicDistArray[FmmBox]{rank==3,rect};

    /**
     * An array of locally essential trees (LETs), one for each place.
     */
    protected val locallyEssentialTrees : DistArray[LocallyEssentialTree](1){rail};

    public static val TIMER_INDEX_MACROSCOPIC : Int = 9;
    /** 
     * A multi-timer for the several segments of a single getEnergy 
     * invocation, indexed by the constants above and in the superclass. 
     */
    public val timer = new Timer(10);

    /** A region representing a cube of 3x3x3 boxes, used for constructing macroscopic multipoles. */
    static val threeCube : Region(3) = (-1..1)*(-1..1)*(-1..1);

    /** A region representing a cube of 9x9x9 boxes, used for interacting with macroscopic multipoles. */
    static val nineCube : Region(3) = (-4..4)*(-4..4)*(-4..4);

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
                    atoms: DistArray[Array[MMAtom](1){rail}](1){rail},
                    numShells : Int) {
        // Periodic FMM always uses ws = 1
        // TODO is it possible to formulate for well-spaced > 1?
        super(density, numTerms, 1, topLeftFront, size, numAtoms, atoms, 0);
        timer.start(TIMER_INDEX_TREE);
        val boxes = constructTree(numLevels, topLevel, numTerms, ws);
        this.boxes = boxes;
        this.lowestLevelBoxes = boxes(numLevels);
        this.numShells = numShells;
        this.locallyEssentialTrees = createLocallyEssentialTrees();
        assignAtomsToBoxes(atoms, boxes(numLevels), offset, lowestLevelDim, size);
        timer.stop(TIMER_INDEX_TREE);
    }
    
    public def calculateEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);
        finish {
            async {
                prefetchPackedAtoms();
            }
            multipoleLowestLevel();
            combineMultipoles();
            combineMacroscopicExpansions();
            transformToLocal();
        }
        val totalEnergy = getDirectEnergy() + getFarFieldEnergy();
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
        val numShells = this.numShells; // TODO shouldn't be necessary XTENLANG-1913
        val numTerms = this.numTerms; // TODO shouldn't be necessary XTENLANG-1913
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val boxes = this.boxes; // TODO shouldn't be necessary XTENLANG-1913
        at (boxes(0).dist(0,0,0)) {
            val macroMultipoles = new Array[MultipoleExpansion](numShells+1);
            val macroLocalTranslations = new Array[LocalExpansion](numShells+1);
            val topLevelBox = boxes(0)(0,0,0);
            macroMultipoles(0) = topLevelBox.multipoleExp;

            var macroTranslation : MultipoleExpansion = new MultipoleExpansion(numTerms);

            // multipoles for shell 1
            for ([i,j,k] in threeCube) {
                val translationVector = Vector3d(i * size,
                                                 j * size,
                                                 k * size);
                val translation = MultipoleExpansion.getOlm(translationVector, numTerms);
                macroTranslation.add(translation);
            }
            macroMultipoles(1) = new MultipoleExpansion(numTerms);
            macroMultipoles(1).translateAndAddMultipole(macroTranslation, macroMultipoles(0));
            //Console.OUT.println("final for 1 = " + macroMultipoles(1));

            // locals for shell 1
            macroLocalTranslations(0) = new LocalExpansion(numTerms);
            for ([i,j,k] in nineCube) {
                if (Math.abs(i) > 1 || Math.abs(j) > 1 || Math.abs(k) > 1) {
                    // inner 27 boxes done at a lower level
                    val translationVector = Vector3d(i * size,
                                                     j * size,
                                                     k * size);
                    val transform = LocalExpansion.getMlm(translationVector, numTerms);
                    macroLocalTranslations(0).add(transform);
                }
            }
            macroLocalTranslations(1) = macroLocalTranslations(0).getMacroscopicParent();

            // remaining shells
            for (var shell: Int = 2; shell <= numShells; shell++) {
                macroTranslation = macroTranslation.getMacroscopicParent();
                macroMultipoles(shell) = new MultipoleExpansion(numTerms);
                macroMultipoles(shell).translateAndAddMultipole(macroTranslation, macroMultipoles(shell-1));
                //Console.OUT.println("final for " + shell + " = " + macroMultipoles(shell));
                macroLocalTranslations(shell) = macroLocalTranslations(shell-1).getMacroscopicParent();
            }

            // now transform and add macroscopic multipoles to local expansion for top level box
            for (var shell: Int = 0; shell <= numShells; shell++) {
                val localExpansion = macroLocalTranslations(shell);
                topLevelBox.localExp.transformAndAddToLocal(localExpansion, macroMultipoles(shell));
            }
            //Console.OUT.println("final for topLevel = " + topLevelBox.localExp);
        }
        timer.stop(TIMER_INDEX_MACROSCOPIC);
    }

    private def assignAtomsToBoxes(atoms: DistArray[Array[MMAtom](1){rail}](1){rail}, lowestLevelBoxes : PeriodicDistArray[FmmBox]{rank==3}, offset : Vector3d, lowestLevelDim : Int, size : Double) {
        //Console.OUT.println("assignAtomsToBoxes");
        val dipole = finish(VectorSumReducer()) {
            ateach (p1 in atoms) {
                var myDipole : Vector3d = Vector3d.NULL;
                val localAtoms = atoms(p1);
                for ([i] in 0..localAtoms.size-1) {
                    val atom = localAtoms(i);
                    val charge = atom.charge;
                    val offsetCentre = atom.centre + offset;
                    myDipole = myDipole + Vector3d(offsetCentre) * charge;
                    val boxIndex = PeriodicFmm3d.getLowestLevelBoxIndex(offsetCentre, lowestLevelDim, size);
                    at(lowestLevelBoxes.dist(boxIndex)) {
                        val remoteAtom = new MMAtom(offsetCentre, charge);
                        val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox;
                        atomic {
                            leafBox.atoms.add(remoteAtom);
                        }
                    }
                }
                offer myDipole;
            }
        };
        this.dipole = dipole;

        // post-prune leaf boxes
        // TODO prune intermediate empty boxes as well
        // TODO pruning before cancel dipole causes NPE on corner 
        //      boxes for small or non-uniform distributions
        finish ateach (boxIndex in lowestLevelBoxes) {
            val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
            if (box.atoms.size() == 0) {
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
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val numTerms = this.numTerms; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach (boxIndex in lowestLevelBoxes) {
            val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox;
            if (leafBox != null) {
                val boxLocation = leafBox.getCentre(size);
                for ([i] in 0..leafBox.atoms.size()-1) {
                    val atom = leafBox.atoms(i);
                    val atomLocation = leafBox.getCentre(size).vector(atom.centre);
                    val atomExpansion = MultipoleExpansion.getOlm(atom.charge, atomLocation, numTerms);
                    leafBox.multipoleExp.add(atomExpansion);
                }
            }
        }

        cancelDipole(dipole);

        timer.stop(TIMER_INDEX_MULTIPOLE);
    }

    def addAtomToLowestLevelBoxAsync(boxIndex : Point(3), offsetCentre : Point3d, charge : Double) {
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val numTerms = this.numTerms; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        async at(lowestLevelBoxes.dist(boxIndex)) {
            val remoteAtom = new MMAtom(offsetCentre, charge);
            val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox;
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
        
        val numLevels = this.numLevels; // TODO shouldn't be necessary XTENLANG-1913
        val topLevel = this.topLevel; // TODO shouldn't be necessary XTENLANG-1913
        val multipoleTransforms = this.multipoleTransforms; // TODO shouldn't be necessary XTENLANG-1913
        val multipoleTranslations = this.multipoleTranslations; // TODO shouldn't be necessary XTENLANG-1913
        val locallyEssentialTrees = this.locallyEssentialTrees; // TODO shouldn't be necessary XTENLANG-1913
        for ([thisLevel] in (topLevel+1)..numLevels) {
            val dim = size / Math.pow2(thisLevel);
            //Console.OUT.println("transform level " + thisLevel);
            val thisLevelBoxes = boxes(thisLevel);
            finish ateach (p1 in Dist.makeUnique()) {
                val myLET = locallyEssentialTrees(p1);
                val combinedVList = myLET.combinedVList;

                val thisLevelMultipoleCopies = myLET.multipoleCopies(thisLevel);
                if (thisLevel == topLevel) {
                    // must fetch top level multipoles synchronously before starting
                    prefetchMultipoles(thisLevel);
                }

                // can fetch next level multipoles asynchronously while computing this level
                val lowerLevel = thisLevel+1;
                if (lowerLevel <= numLevels) {
                    async prefetchMultipoles(lowerLevel);
                }

                for ([x1,y1,z1] in thisLevelBoxes.dist(here)) async {
                    val box1 = thisLevelBoxes(x1,y1,z1);
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for (p in vList) {
                            val boxIndex = vList(p);
                            // force on the multipole value
                            val box2MultipoleExp = thisLevelMultipoleCopies(boxIndex);
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(boxIndex);
                                val translateP = Point.make([here.id, thisLevel, translation(0), translation(1), translation(2)]);
                                val transform21 = multipoleTransforms(translateP);
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                        val box1Parent = box1.parent;
                        val box1ParentExp = at (box1Parent) {box1Parent().localExp};
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, box1.x%2, box1.y%2, box1.z%2]));
                        box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_TRANSFORM);
    }

    /**
     * Fetch all multipole expansions required for transform to local at this place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchMultipoles(level : Int) : void {
        val myLET = locallyEssentialTrees(here.id);
        val combinedVList = myLET.combinedVList(level);
        val thisLevelCopies = myLET.multipoleCopies(level);
        val thisLevelBoxes = boxes(level);

        val vListPlaces = new HashMap[Int,ArrayList[Point(3)]](26); // a place may have up to 26 immediate neighbours
        
        // separate the vList into partial lists stored at each nearby place
        for (p in combinedVList) {
            val boxIndex = combinedVList(p);
            val placeId = thisLevelBoxes.periodicDist(boxIndex).id;
            var vListForPlace : ArrayList[Point(3)] = vListPlaces.getOrElse(placeId, null);
            if (vListForPlace == null) {
                vListForPlace = new ArrayList[Point(3)]();
                vListPlaces.put(placeId, vListForPlace);
            }
            vListForPlace.add(boxIndex);
        }

        // retrieve the partial list for each place and store into my LET
        finish for(placeEntry in vListPlaces.entries()) async {
            val placeId = placeEntry.getKey();
            val vListForPlace = placeEntry.getValue();
            val vListArray = vListForPlace.toArray();
            val multipolesForPlace = at (Place.place(placeId)) { getMultipolesForBoxList(thisLevelBoxes, vListArray)};
            for ([i] in 0..vListArray.size-1) {
                thisLevelCopies(vListArray(i)) = multipolesForPlace(i);
            }
        }
    }

    /**
     * Given a level and a list of box indexes (as Point(3)) stored
     * at a single place, returns an Array, each element of which 
     * is in turn a MultipoleExpansion for the box.
     */
    private static def getMultipolesForBoxList(thisLevelBoxes : PeriodicDistArray[FmmBox]{rank==3}, boxList : Array[Point(3)](1){rail}) {
        val multipoleList = new Array[MultipoleExpansion](boxList.size, 
                                                            (i : Int) => 
                                                                getMultipoleExpansionLocalCopy(thisLevelBoxes,
                                                                                     boxList(i)(0), 
                                                                                     boxList(i)(1),
                                                                                     boxList(i)(2))
                                                            );
        return multipoleList;
    }

    /**
     * Fetch all atoms required for direct calculations at each place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchPackedAtoms() : void {
        timer.start(TIMER_INDEX_PREFETCH);
        val locallyEssentialTrees = this.locallyEssentialTrees; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach (p1 in locallyEssentialTrees) {
            val myLET = locallyEssentialTrees(p1);
            val myCombinedUList = myLET.combinedUList;
            val packedAtoms = myLET.packedAtoms;

            val uListPlaces = new HashMap[Int,ArrayList[Point(3)]](26); // a place may have up to 26 immediate neighbours
            
            // separate the uList into partial lists stored at each nearby place
            for (p in myCombinedUList) {
                val boxIndex = myCombinedUList(p);
                val placeId = lowestLevelBoxes.periodicDist(boxIndex).id;
                var uListForPlace : ArrayList[Point(3)] = uListPlaces.getOrElse(placeId, null);
                if (uListForPlace == null) {
                    uListForPlace = new ArrayList[Point(3)]();
                    uListPlaces.put(placeId, uListForPlace);
                }
                uListForPlace.add(boxIndex);
            }

            // retrieve the partial list for each place and store into my LET
            finish for(placeEntry in uListPlaces.entries()) async {
                val placeId = placeEntry.getKey();
                val uListForPlace = placeEntry.getValue();
                val uListArray = uListForPlace.toArray();
                val packedForPlace = at (Place.place(placeId)) { getPackedAtomsForBoxList(uListArray)};
                for ([i] in 0..uListArray.size-1) {
                    myLET.packedAtoms(uListArray(i)) = packedForPlace(i);
                }
            }
        }
        timer.stop(TIMER_INDEX_PREFETCH);
        //Console.OUT.println("done prefetch");
    }

    /**
     * Given a list of box indexes (as Point(3) stored at a single
     * place, returns an Array, each element of which is in turn
     * an Array of MMAtom.PackedRepresentation containing the 
     * packed atoms for each box.
     */
    private def getPackedAtomsForBoxList(boxList : Array[Point(3)](1){rail}) {
        val packedAtomList = new Array[Array[MMAtom.PackedRepresentation](1){rail}](boxList.size, 
                                                            (i : Int) => 
                                                                getPackedAtomsForBox(boxList(i)(0), 
                                                                                     boxList(i)(1),
                                                                                     boxList(i)(2))
                                                            );
        return packedAtomList;
    }

    private def constructTree(numLevels : Int, topLevel : Int, numTerms : Int, ws : Int) 
      : Array[PeriodicDistArray[FmmBox]{rank==3,rect}](1){rail} {
        val boxArray = new Array[PeriodicDistArray[FmmBox]{rank==3,rect}](numLevels+1);
        for ([thisLevel] in topLevel..numLevels) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(3){rect} = (0..levelDim-1) * (0..levelDim-1) * (0..levelDim-1);
            val thisLevelDist = MortonWorldDist.make(thisLevelRegion);
            boxArray(thisLevel) = PeriodicDistArray.make[FmmBox](thisLevelDist);
            //Console.OUT.println("level " + thisLevel + " dist: " + thisLevelDist);
        }

        for ([thisLevel] in topLevel..numLevels-1) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelBoxes = boxArray(thisLevel);
            finish ateach ([x,y,z] in thisLevelBoxes) {
                val box = new FmmBox(thisLevel, x, y, z, numTerms, PeriodicFmm3d.getParentForChild(boxArray, thisLevel, topLevel, x,y,z));
                createVList(box, ws);
                thisLevelBoxes(x,y,z) = box;
            }
        }

        val lowestLevelBoxes = boxArray(numLevels);
        finish ateach ([x,y,z] in lowestLevelBoxes) {
            val box = new FmmLeafBox(numLevels, x, y, z, numTerms, PeriodicFmm3d.getParentForChild(boxArray, numLevels, topLevel, x,y,z));
            createUList(box, ws);
            createVList(box, ws);
            lowestLevelBoxes(x,y,z) = box;
        }

        return boxArray;
    }

    /**
     * Creates the U-list of <code>box</code>.
     * The U-list consists of all leaf boxes not well-separated from <code>box</code>.
     */
    private static def createUList(box : FmmLeafBox, ws : Int) {
        // interact with "left half" of uList i.e. only boxes with x<=box.x
        val uList = new ArrayList[Point(3)]();
        for ([x] in box.x-ws..box.x) {
            for ([y] in box.y-ws..box.y+ws) {
                for ([z] in box.z-ws..box.z+ws) {
                    if (x < box.x || (x == box.x && y < box.y) || (x == box.x && y == box.y && z < box.z)) {
                        uList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setUList(uList.toArray());
    }

    /**
     * Creates the V-list of <code>box</code>.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent of <code>box</code>.
     */
    private static def createVList(box : FmmBox, ws : Int) {
        val xOffset = box.x%2 == 1 ? -1 : 0;
        val yOffset = box.y%2 == 1 ? -1 : 0;
        val zOffset = box.z%2 == 1 ? -1 : 0;
        val vList = new ArrayList[Point(3)]();
        for ([x] in box.x-2*ws+xOffset..box.x+2*ws+1+xOffset) {
            for ([y] in box.y-2*ws+yOffset..box.y+2*ws+1+yOffset) {
                for ([z] in box.z-2*ws+zOffset..box.z+2*ws+1+zOffset) {
                    if (box.wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setVList(vList.toArray());
    }

    /**
     * Gets sum of direct (pairwise) energy for all pairs of atoms
     * in non-well-separated boxes. This operations requires only
     * that atoms have already been assigned to boxes, and so can 
     * be done in parallel with other steps of the algorithm.
     */
    def getDirectEnergy() : Double {
        //Console.OUT.println("direct");
        timer.start(TIMER_INDEX_DIRECT);

        val locallyEssentialTrees = this.locallyEssentialTrees; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        val directEnergy = finish (SumReducer()) {
            ateach (p1 in locallyEssentialTrees) {
                val myLET = locallyEssentialTrees(p1);
                val packedAtoms = myLET.packedAtoms;
                var thisPlaceEnergy : Double = 0.0;
                for ([x1,y1,z1] in lowestLevelBoxes.dist(here)) {
                    val box1 = lowestLevelBoxes(x1,y1,z1) as FmmLeafBox;
                    if (box1 != null) {
                        for ([atomIndex1] in 0..box1.atoms.size()-1) {
                            val atom1 = box1.atoms(atomIndex1);

                            // direct calculation with all atoms in same box
                            for ([sameBoxAtomIndex] in 0..atomIndex1-1) {
                                val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                                val pairEnergy = atom1.charge * sameBoxAtom.charge / atom1.centre.distance(sameBoxAtom.centre);
                                thisPlaceEnergy += pairEnergy;
                            }

                            // direct calculation with all atoms in non-well-separated boxes
                            val uList = box1.getUList();
                            for (p in uList) {
                                val boxIndex2 = uList(p);
                                val boxAtoms = packedAtoms(boxIndex2);
                                if (boxAtoms != null) {
                                    for ([otherBoxAtomIndex] in 0..(boxAtoms.size-1)) {
                                        val atom2Packed = boxAtoms(otherBoxAtomIndex);
                                        thisPlaceEnergy += atom1.charge * atom2Packed.charge / atom1.centre.distance(atom2Packed.centre);
                                    }
                                }
                            }
                        }
                    }
                }
                offer thisPlaceEnergy;
            }
        };
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

    private static def getParentForChild(boxes : Array[PeriodicDistArray[FmmBox]{rank==3,rect}](1){rail}, level : Int, topLevel : Int, x : Int, y : Int, z : Int) : GlobalRef[FmmBox] {
        if (level == topLevel)
            return GlobalRef[FmmBox](null);
        return (at (boxes(level-1).dist(x/2, y/2, z/2)) {GlobalRef[FmmBox](boxes(level-1)(x/2, y/2, z/2))});
    }

    /**
     * Creates the locally essential tree at each place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private def createLocallyEssentialTrees() : DistArray[LocallyEssentialTree]{rank==1} {
        val locallyEssentialTrees = DistArray.make[LocallyEssentialTree](Dist.makeUnique(), (Point)=> null);
        val topLevel = this.topLevel; // TODO shouldn't be necessary XTENLANG-1913
        val numLevels = this.numLevels; // TODO shouldn't be necessary XTENLANG-1913
        val boxes = this.boxes; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach ([p1] in locallyEssentialTrees) {
            val uMin = new Array[Int](3, Int.MAX_VALUE);
            val uMax = new Array[Int](3, Int.MIN_VALUE);
            val combinedUSet = new HashSet[Point(3)]();
            for ([x,y,z] in lowestLevelBoxes.dist(here)) {
                val box1 = lowestLevelBoxes(x,y,z) as FmmLeafBox;
                if (box1 != null) {
                    val uList = box1.getUList();
                    for (p in uList) {
                        val boxIndex2 = uList(p);
                        if (combinedUSet.add(boxIndex2)) {
                            for ([i] in 0..2) {
                                uMin(i) = Math.min(uMin(i), boxIndex2(i));
                                uMax(i) = Math.max(uMax(i), boxIndex2(i));        
                            }
                        }
                    }
                }
            }
            val combinedUList = new Array[Point(3)](combinedUSet.size());
            var j : Int = 0;
            for (boxIndex in combinedUSet) {
                combinedUList(j++) = boxIndex;
            }

            val combinedVList = new Array[Array[Point(3)](1){rail}](numLevels+1);
            val vListMin = new Array[Array[Int](1){rail}](numLevels+1);
            val vListMax = new Array[Array[Int](1){rail}](numLevels+1);
            for ([thisLevel] in topLevel..numLevels) {
                val vMin = new Array[Int](3, Int.MAX_VALUE);
                val vMax = new Array[Int](3, Int.MIN_VALUE);
                //Console.OUT.println("create combined V-list for level " + thisLevel + " at " + here);
                val combinedVSet = new HashSet[Point(3)]();
                val thisLevelBoxes = boxes(thisLevel);
                for ([x,y,z] in thisLevelBoxes.dist(here)) {
                    val box1 = thisLevelBoxes(x,y,z);
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for (p in vList) {
                            val boxIndex2 = vList(p);
                            if (combinedVSet.add(boxIndex2)) {
                                for ([i] in 0..2) {
                                    vMin(i) = Math.min(vMin(i), boxIndex2(i));
                                    vMax(i) = Math.max(vMax(i), boxIndex2(i));        
                                }
                            }
                        }
                    }
                }
                //Console.OUT.println("done " + combinedVSet.size());
                vListMin(thisLevel) = vMin;
                vListMax(thisLevel) = vMax;
                combinedVList(thisLevel) = new Array[Point(3)](combinedVSet.size());
                var i : Int = 0;
                for (boxIndex in combinedVSet) {
                    combinedVList(thisLevel)(i++) = boxIndex;
                }
            }
            locallyEssentialTrees(p1) = new LocallyEssentialTree(combinedUList,
                                                                 combinedVList,
                                                                 uMin,
                                                                 uMax,
                                                                 vListMin,
                                                                 vListMax);
        }
        return locallyEssentialTrees;
    }

    /**
     * @return a local copy at the current place of a box's multipole expansion
     */
    /**
     * @return a local copy at the current place of a box's multipole expansion
     */
    private static def getMultipoleExpansionLocalCopy(thisLevelBoxes : PeriodicDistArray[FmmBox]{rank==3}, x : Int, y : Int, z : Int) : MultipoleExpansion {
        return at (thisLevelBoxes.periodicDist(x,y,z)) {thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null};
    }

    private def getPackedAtomsForBox(x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox;
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
        for (var level: Int = numLevels; level > topLevel; level--) {
            val thisLevel = level;
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelBoxes = boxes(thisLevel);
            val multipoleTranslations = this.multipoleTranslations; // TODO shouldn't be necessary XTENLANG-1913
            finish ateach (boxIndex in thisLevelBoxes) {
                val child = thisLevelBoxes(boxIndex);
                if (child != null) {
                    val childExp = child.multipoleExp;
                    val parent = child.parent;
                    at (parent) {
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, (child.x+1)%2, (child.y+1)%2, (child.z+1)%2]));
                        parent().multipoleExp.translateAndAddMultipole(shift, childExp);
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
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelBoxes = this.lowestLevelBoxes; // TODO shouldn't be necessary XTENLANG-1913
        val farFieldEnergy = finish(SumReducer()) {
            ateach (boxIndex1 in lowestLevelBoxes) {
                val box1 = lowestLevelBoxes(boxIndex1) as FmmLeafBox;
                if (box1 != null) {
                    var thisBoxEnergy : Double = 0.0;
                    for ([atomIndex1] in 0..box1.atoms.size()-1) {
                        val atom1 = box1.atoms(atomIndex1);
                        val box1Centre = atom1.centre.vector(box1.getCentre(size));
                        val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                        thisBoxEnergy += farFieldEnergy;
                    }
                    offer thisBoxEnergy;
                }
            }
        };
        timer.stop(TIMER_INDEX_FARFIELD);

        return farFieldEnergy / 2.0;
    }

    private static def getLowestLevelBoxIndex(offsetCentre : Point3d, lowestLevelDim : Int, size : Double) : Point(3) {
        return  Point.make((offsetCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as Int);
    }

    static struct VectorSumReducer implements Reducible[Vector3d] {
        public def zero() = Vector3d.NULL;
        public def apply(a:Vector3d, b:Vector3d) = (a + b);
    }

}

