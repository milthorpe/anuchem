package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed Particle Mesh Ewald implementation.
 * @author milthorpe
 */
public class TestPME extends TestElectrostatic {
    public global def sizeOfCentralCluster() : Double = 2.0;

    public def this(numAtoms : Int) {
        super(numAtoms);
    }

    public static def main(args : Rail[String]!) {
        var numAtoms : Int;
        var ewaldCoefficient : Double = 0.35;
        var cutoff : Double = 10.0;
        var gridSize : Int = 64;
        var splineOrder : Int = 4;
        if (args.length > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.length > 1) {
                ewaldCoefficient = Double.parseDouble(args(1));
                if (args.length > 2) {
                    cutoff = Double.parseDouble(args(2));
                    if (args.length > 3) {
                        gridSize = Int.parseInt(args(3));
                        if (args.length > 4) {
                            splineOrder = Int.parseInt(args(4));
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPME numAtoms [ewaldCoefficient] [cutoff] [gridSize] [splineOrder]");
            return;
        }

        if (splineOrder > gridSize) {
            Console.ERR.println("TestPME: splineOrder must not be greater than gridSize");
            return;
        }
        new TestPME(numAtoms).test(ewaldCoefficient, cutoff, gridSize, splineOrder);
    }

    public def test(ewaldCoefficient : Double, cutoff : Double, gridSize : Int, splineOrder : Int) {

        val edges = [Vector3d(SIZE, 0.0, 0.0), Vector3d(0.0, SIZE, 0.0), Vector3d(0.0, 0.0, SIZE)];
        val g = gridSize;
        val gridSizes = ValRail.make[Int](3, (Int) => g);

        Console.OUT.println("Testing PME for " + numAtoms + " particles."
            + "\nBox edges: " + edges
            + "\nGrid size: " + gridSize
            + "\nspline order: " + splineOrder + " Beta: " + ewaldCoefficient + " Cutoff: " + cutoff);

        val atoms = generateAtoms();
        val pme = new PME(edges, gridSizes, atoms, splineOrder, ewaldCoefficient, cutoff);
        val energy = pme.getEnergy();
        Console.OUT.println("energy = " + energy);

        logTime("Divide",            PME.TIMER_INDEX_DIVIDE,        pme.timer);
        logTime("Direct",            PME.TIMER_INDEX_DIRECT,        pme.timer);
        logTime("Self energy",       PME.TIMER_INDEX_SELF,          pme.timer);
        logTime("Grid charges",      PME.TIMER_INDEX_GRIDCHARGES,   pme.timer);
        logTime("Inverse FFT",       PME.TIMER_INDEX_INVFFT,        pme.timer);
        logTime("ThetaRecConvQ",     PME.TIMER_INDEX_THETARECCONVQ, pme.timer);
        logTime("Reciprocal energy", PME.TIMER_INDEX_RECIPROCAL,    pme.timer);
        logTime("Total",             PME.TIMER_INDEX_TOTAL,         pme.timer);

        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("cf. Direct calculation", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        // direct error comparison is only useful if there is a huge empty border around the particles
        val error = directEnergy - energy;
        Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));
    }
}

