/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.pme;

import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.util.ArrayList;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.GromacsStructureFileReader;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed Particle Mesh Ewald implementation using
 * the "bigwater" box - an 80 Angstrom cube filled with SPC water.
 * atom coordinates are read from a GROMACS structure file.
 * @see H. Berendsen, J. Grigera, and T. Straatsma,
 *      "The missing term in effective pair potentials,"
 *      J. Phys. Chem., vol. 91, pp. 6269-71, Nov 1987.
 * @author milthorpe
 */
public class TestPMEWaterBox extends TestElectrostatic {
    public static ITERS = 10;
    var size:Double = 80.0;
    public def sizeOfCentralCluster() { return size; }
    public def boxSize() { return size; }

    public static def main(args:Rail[String]) {
        var structureFileName : String;
        var ewaldCoefficient : Double = 0.35;
        var cutoff : Double = 10.0;
        var gridSize : Long = 72;
        var splineOrder : Int = 4n;
        if (args.size > 0) {
            structureFileName = args(0);
            if (args.size > 1) {
                ewaldCoefficient = Double.parseDouble(args(1));
                if (args.size > 2) {
                    cutoff = Double.parseDouble(args(2));
                    if (args.size > 3) {
                        gridSize = Long.parseLong(args(3));
                        if (args.size > 4) {
                            splineOrder = Int.parseInt(args(4));
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: pme myStructureFile.gro [ewaldCoefficient] [cutoff] [gridSize] [splineOrder]");
            return;
        }

        if (splineOrder > gridSize) {
            Console.ERR.println("pme: splineOrder must not be greater than gridSize");
            return;
        }
        new TestPMEWaterBox().test(structureFileName, ewaldCoefficient, cutoff, gridSize, splineOrder);
    }

    public def test(structureFileName : String, ewaldCoefficient : Double, cutoff : Double, gridSize : Long, splineOrder : Int) {

        val gmxFileReader = new GromacsStructureFileReader(structureFileName);
        val molecule = gmxFileReader.readMolecule();
        this.size = gmxFileReader.boxEdges.i; // TODO assume cubic box

        val edges = [Vector3d(size, 0.0, 0.0), Vector3d(0.0, size, 0.0), Vector3d(0.0, 0.0, size)];
        val g = gridSize;
        val gridSizes = new Rail[Long](3, g);

        Console.OUT.println("Testing PME with structure file " + structureFileName
            + "\nBox edges: " + edges(0) + "," + edges(1) + "," + edges(2)
            + "\nGrid size: " + gridSize
            + "\nspline order: " + splineOrder + " Beta: " + ewaldCoefficient + " Cutoff: " + cutoff);

        val molAtoms = molecule.getAtoms();
        Console.OUT.println("read " + molAtoms.size() + " atoms.");

        val tempAtoms = DistArray.make[ArrayList[MMAtom]](Dist.makeUnique(), (Point) => new ArrayList[MMAtom]());
        finish for (var i:Long = 0; i < molAtoms.size(); i++) {
            val atom = molAtoms(i);
            val centre = atom.centre;

            // normalise coordinates to fit inside box
            val x = -size/2.0 + ((centre.i < 0.0) ? centre.i + size : (centre.i >= size) ? centre.i - size : centre.i);
            val y = -size/2.0 + ((centre.j < 0.0) ? centre.j + size : (centre.j >= size) ? centre.j - size : centre.j);
            val z = -size/2.0 + ((centre.k < 0.0) ? centre.k + size : (centre.k >= size) ? centre.k - size : centre.k);

            val charge = atom.charge;
            val p = getPlaceId(x, y, z);
            if (p >= 0 && p < Place.MAX_PLACES) {
                at(Place.place(p)) async {
                    atom.centre = Point3d(x,y,z);
                    //Console.OUT.println(atom);
                    atomic { tempAtoms(p).add(atom); }
                }
            } else {
                Console.ERR.println("could not map atom to place: " + atom.centre);
            }
        }
        val atoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique(), ([p] : Point) => tempAtoms(p).toRail());

        val pme = new PME_SPMD(edges, gridSizes, atoms, splineOrder, ewaldCoefficient, cutoff);
        pme.setup();
        val energy = pme.getEnergy();
        for (i in 1..ITERS) {
            pme.getEnergy();
        }
        Console.OUT.println("energy = " + energy);

        logTime("Grid charges",      PME.TIMER_INDEX_GRIDCHARGES,   pme.timer);
        logTime("Inverse FFT",       PME.TIMER_INDEX_INVFFT,        pme.timer);
        logTime("ThetaRecConvQ",     PME.TIMER_INDEX_THETARECCONVQ, pme.timer);
        logTime("Reciprocal energy", PME.TIMER_INDEX_RECIPROCAL,    pme.timer);
        logTime("Self energy",       PME.TIMER_INDEX_SELF,          pme.timer);
        logTime("Direct",            PME.TIMER_INDEX_DIRECT,        pme.timer);
        logTime("Total",             PME.TIMER_INDEX_TOTAL,         pme.timer);
        logTime("Prefetch",          PME.TIMER_INDEX_PREFETCH,      pme.timer);
        logTime("Setup",             PME.TIMER_INDEX_SETUP,         pme.timer);
    }
}

