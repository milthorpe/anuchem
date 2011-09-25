/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */

import x10.compiler.NonEscaping;
import x10.compiler.Inline;
import au.edu.anu.util.Timer;

/**
 * A simple 2D lattice Boltzmann simulation of flow in a lid-driven square cavity.
 * Based on Fortran/GA example distributed with the Global Arrays library.
 * @see http://www.emsl.pnl.gov/docs/global/
 * @see Palmer, B.J. and Rector, D.R. (2000). "Lattice Boltzmann algorithm for
 *   simulating thermal flow in compressible fluids". J. Comp. Phys. 161 1-20
 * @see Palmer, B.J. and Nieplocha, J. (2002). "Efficient algorithms for ghost 
 *   cell updates on two classes of MPP architectures". Proc. PDCS 2002
 *   http://www.emsl.pnl.gov/docs/global/papers/ghosts.pdf
 * @author milthorpe 09/2011
 * 23/09 started 20:45
 */
public class Boltzmann {
    public static NSIZE = 129;
    public static GHOST_WIDTH = 1;

    public static NSTEPS = 30;

    public static VISCOSITY = 1.0;

    /** time increment*/
    public static DELTA_T = 1.0;

    /** the dimensions of the square cavity */
    public static XMAX = 256.0;

    /** temperature of fluid */
    public static TMPRTR0 = 1.0;

    /** uniform density of fluid in cavity */
    public static RHO0 = 2.7;

    /** velocity of lid */
    public static UXBC = 0.5;

    /** the ideal gas constant (eq. 4.1, 4.2, 4.9, 4.10) */
    public static RGAS = 1.0;

    /** van der Waals parameter a (eq. 4.9) */
    public static A_VDW = 0.0;

    /** van der Waals parameter b (eq. 4.10) */
    public static B_VDW = 0.0;

    /** The region of the full lattice, including boundary points. */
    val latticeRegion = 0..(NSIZE-1) * 0..(NSIZE-1);

    /** The distributed of the full lattice. */
    val latticeDist = Dist.makeBlockBlock(latticeRegion, 0, 1);

    /** The region of interior points of the lattice. */
    val interiorRegion = 1..(NSIZE-2) * 1..(NSIZE-2);

    /** The lattice Boltzmann distribution functions at each site. maps to g_fg*/
    val lb:DistArray[LBDist](3);

    /** The properties (density, momentum etc.) at each site. maps to g_fld*/
    val lbProps:DistArray[LBProperties](2);

    /** The boundary condition data. maps to g_bc*/
    val boundary:DistArray[Int](2);

    /** Total mass */
    var rtot:Double;

    /** Total x-momentum */
    var uxtot:Double;

    /** Total y-momentum */
    var uytot:Double;

    /** Speed - magnitude of velocity for e_i */
    val cspd:Double;

    public val timer = new Timer(2); // maps to tstats

    val ffb = new Array[Double](9);
    val ffc = new Array[Double](9);
    val ffd = new Array[Double](9);

    /** The velocity vectors e_i */
    static ei = initEi();

    static hash = initHash();
    static ihash = initIHash();

    public def this() {
        val ux0 = 0.0;
        val uy0 = 0.0;

        val delta_x = XMAX/((NSIZE-1) as Double);
        val cspd = Math.sqrt(2.0)*delta_x/DELTA_T;
        this.cspd = cspd;
        val tau_rho = 6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5;
        Console.OUT.println("Value of tau_rho is " + tau_rho);

        //Console.OUT.println("latticeDist = " + latticeDist);

        // initialize boundary conditions
        val boundaryCond = ([i,j]:Point(2)) => { (i==1 || i==NSIZE || j==1) ? 1 : (j==NSIZE) ? 2 : 0 };
        val boundary = DistArray.make[Int](latticeDist, boundaryCond, GHOST_WIDTH);
        boundary.updateGhosts();
        this.boundary = boundary;

        // initialize distribution of density and velocities
        val lbProps = DistArray.make[LBProperties](latticeDist, 
            (Point) => new LBProperties(
                RHO0, // + 0.0*cos(2.0d00*pi*dble(j-1) / dble(size(2)-1))                                    
                ux0, uy0,
                RHO0 * RGAS * TMPRTR0 / (1.0 - B_VDW*RHO0) - A_VDW*RHO0*RHO0,
                6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5),
            GHOST_WIDTH);
        this.lbProps = lbProps;
        rtot = RHO0*lbProps.region.size();

        initializeLatticeParameters();

        val lbDist = Dist.makeBlockBlock(latticeRegion * 0..8, 0, 1);
        //Console.OUT.println("lbDist = " + lbDist);
        val lb = DistArray.make[LBDist](lbDist, (Point) => new LBDist(), GHOST_WIDTH);
        this.lb = lb;
    }
    

    public def boltzmann() {
        // initialize equilibrium and actual distributions
        equilibrium();
        finish ateach([ii,jj] in latticeDist | interiorRegion) {
            for (i in 0..8) {
                lb(ii,jj,i).actual = lb(ii,jj,i).equilibrium;
            }
        }
        Console.OUT.println("Total mass is " + rtot);
        Console.OUT.println("Total x-momentum is " + uxtot);
        Console.OUT.println("Total y-momentum is " + uytot);

        for (i in 1..NSTEPS) {
            timestep();
            if (i % 200 == 0) {
                Console.OUT.println("Completed step " + i);
                Console.OUT.println("Total density is " + rtot);
                Console.OUT.println("Total x-momentum is " + uxtot);
                Console.OUT.println("Total y-momentum is " + uytot);
                vorticity();
                printData();
            }
        }
        vorticity();
        printData();
        Console.OUT.println("Time in lattice Boltzmann updates :" + (timer.total(0) as Double) / 1e9);
        Console.OUT.println("Time in ghost cell updates        :" + (timer.total(1) as Double) / 1e9);
    }

    private @NonEscaping def initializeLatticeParameters() {
        val rdim = 4.0;
        val c2 = cspd*cspd;
        val b = 24.0;

        // TODO make ffb, ffc, ffd PlaceLocalHandles
        ffb(0) = 0.0;
        ffc(0) = 0.0;
        ffd(0) = -1.0 / c2 - 4.0 * rdim / (2.0 * b * c2);
        for (i in 1..4) {
          ffb(i) = 4.0 *rdim/(b*c2);
          ffc(i) = 4.0 *rdim*(rdim+2.0)/(2.0*b*c2*c2);
          ffd(i) = -4.0 *rdim/(2.0*b*c2);
        }
        for (i in 5..8) {
          ffb(i) = rdim/(b*c2);
          ffc(i) = rdim*(rdim+2.0)/(2.0*b*c2*c2);
          ffd(i) = -rdim/(2.0*b*c2);
        }
    }

    private static def initHash(): Array[Int](2){rect} {
        val hash = new Array[Int](-1..1 * -1..1);
        hash(0,0)   = 1;
        hash(1,0)   = 2;
        hash(-1,0)  = 3;
        hash(0,1)   = 4;
        hash(0,-1)  = 5;
        hash(1,1)   = 6;
        hash(-1,-1) = 7;
        hash(1,-1)  = 8;
        hash(-1,1)  = 9;
        return hash;
    }

    private static def initIHash(): Array[Int](2){rect} {
        val ihash = new Array[Int](-1..1 * -1..1);
        ihash(0,0)   = 1;
        ihash(1,0)   = 3;
        ihash(-1,0)  = 2;
        ihash(0,1)   = 5;
        ihash(0,-1)  = 4;
        ihash(1,1)   = 7;
        ihash(-1,-1) = 6;
        ihash(1,-1)  = 9;
        ihash(-1,1)  = 8;
        return ihash;
    }

    private static def initEi(): Array[Int](2){rect,zeroBased} {
        val ei = new Array[Int](0..8 * 0..1);
        ei(0,0) =  1; ei(0,1) =  0;
        ei(1,0) = -1; ei(1,1) =  0;
        ei(2,0) =  0; ei(2,1) =  1;
        ei(3,0) =  0; ei(3,1) = -1;
        ei(4,0) =  1; ei(4,1) =  1;
        ei(5,0) = -1; ei(5,1) = -1;
        ei(6,0) =  1; ei(6,1) = -1;
        ei(7,0) = -1; ei(7,1) =  1;
        ei(8,0) =  0; ei(8,1) =  0;
        return ei;
    }

    /** Update values of density, momentum and pressure */
    private def properties() {
        val cspd2 = cspd / Math.sqrt(2.0);
		val rtot = new Accumulator[Double](SumReducer());
		val uxtot = new Accumulator[Double](SumReducer());
		val uytot = new Accumulator[Double](SumReducer());
        finish ateach (p in Dist.makeUnique()) {
            val lbInterior = (latticeDist.region && interiorRegion) as Region(2){rect};
            for([ii,jj] in lbInterior) {
                val props = lbProps(ii,jj);

                // evaluate density and momentum
                for (i in 0..8) {
                    val j = i-1;
                    val ex:Double;
                    val ey:Double;
                    if (j==-1) {
                        ex = 0.0;
                        ey = 0.0;
                    } else {
                        ex = cspd2*ei(j,0);
                        ey = cspd2*ei(j,1);
                    }
                    props.density += lb(ii,jj,i).actual;
                    props.px += ex * lb(ii,jj,i).actual;
                    props.py += ey * lb(ii,jj,i).actual;
                }

                // note: momentum accumulated *before* dividing by density
                uxtot <- props.px;
                uytot <- props.py;

                props.px /= props.density;
                props.py /= props.density;

    /*
                // TODO boundary modifications are commented out in GA code
                if (boundary(ii,jj) == 1) {
                    // sides
                    props.px = 0.0;
                    props.py = 0.0;
                } else if (boundary(ii,jj) == 2) {
                    // lid
                    props.density = RHO0;
                    props.px = UXBC;
                    props.py = 0.0;
                }
    */

                val rho = props.density;
                // evaluate pressure (eq. 4.9)
                props.pressure = rho*RGAS*TMPRTR0/(1.0-B_VDW*rho) - A_VDW*rho*rho;

                rtot <- rho;
            }
        }
        this.rtot = rtot();
        this.uxtot = uxtot();
        this.uytot = uytot();
    }

    /** Update equilibrium distributions */
    private def equilibrium() {
        val rdim = 4.0;
        val b = 24.0;
        val c2 = cspd*cspd;
        val cspdi = Math.sqrt(2.0)/cspd;
        val cspd2 = 1.0/cspdi;

        finish ateach (p in Dist.makeUnique()) {
            val ffa = new Array[Double](9);
            val lbInterior = (latticeDist.region && interiorRegion) as Region(2){rect};
            for([ii,jj] in lbInterior) {
                for (j in 0..8) {
                    val ex:Double;
                    val ey:Double;
                    if (j==0) {
                        ex = 0.0;
                        ey = 0.0;
                    } else {
                        ex = cspd2*ei(j-1,0);
                        ey = cspd2*ei(j-1,1);
                    }
                    val props = lbProps(ii,jj);
                    val aa = (rdim/(b*c2)) * props.pressure;
                    val aa0 = props.density - b*aa;
                    if (j == 1) {
                        ffa(j) = aa0 + 4.0*aa;
                    } else if (j < 5) {
                        ffa(j) = 4.0*aa;
                    } else {
                        ffa(j) = aa;
                    }

                    lb(ii,jj,j).equilibrium = ffa(j) 
                       + props.density * ffb(j)*(props.px*ex + props.py*ey)
                       + ffc(j) * 
                            (props.px*props.px * ex*ex
                              + 2.0*props.px*props.py*ex*ey
                              + props.py*props.py * ey*ey)
                       + ffd(j) * (props.px*props.px + props.py*props.py);
                }
            }
        }
    }

    private def timestep() {
        timer.start(0);

        val lb = this.lb;
        val latticeDist = this.latticeDist;
        val interiorRegion = this.interiorRegion;

        // make backup copies of distribution and update ghost cells
        //finish ateach ([ii,jj] in latticeDist | interiorRegion) {
        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbInterior = (lb.dist(here) && (interiorRegion * 0..8)) as Region(3){rect};
            for ([ii,jj,i] in lbInterior) {
                lb(ii,jj,i).temp = lb(ii,jj,i).actual;
            }

        }

        timer.start(1);
        lb.updateGhosts();
        timer.stop(1);

        // Perform streaming operation
        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbInterior = (latticeDist.region && interiorRegion) as Region(2){rect};
            val patch = new Array[Double](-1..1 * -1..1 * 0..8);
            for ([ii,jj] in lbInterior) {
                if (boundary(ii,jj) == 0) {
                    for (i in 1..8) {
                        val k = i-1;
                        val ix = ei(k,0); //Math.round(ei(k,0)) as Int;
                        val iy = ei(k,1); //Math.round(ei(k,1)) as Int;
                        lb(ii,jj,i).actual = lb(ii-ix,jj-iy,i).temp;
                    }
                } else {
                    getPatch(ii, jj, patch); // maps to fgp
                    for (i in 1..8) {
                        val k = i-1;
                        val ix = ei(k,0); //Math.round(ei(k,0)) as Int;
                        val iy = ei(k,1); //Math.round(ei(k,1)) as Int;
                        lb(ii,jj,i).actual = patch(-ix,-iy,i);
                    }   
                }
            }
        }

        properties();

        // perform relaxation
        equilibrium();

        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbInterior = (lb.dist(here) && (interiorRegion * 0..8)) as Region(3){rect};
            for ([ii,jj,i] in lbInterior) {  
                if (lbProps(ii,jj).t_rho > 0.0) {
                    lb(ii,jj,i).actual = lb(ii,jj,i).temp
                             - (lb(ii,jj,i).temp - lb(ii,jj,i).equilibrium) / lbProps(ii,jj).t_rho;
                }
                // TODO boundary constant condition is commented out in GA code
    //            if (bc(ii,jj) == 2) {
    //                lb(ii,jj,i).actual = lb(ii,jj,i).equilibrium;
    //            }
            }
        }
        timer.stop(0);
    }

    private @Inline def getPatch(ii:Int, jj:Int, patch:Array[Double](3){rect}) {
        // TODO
    }

    private def vorticity() {

    }

    private def printData() {

    }

	static struct SumReducer implements Reducible[Double] {
        public def zero():Double = 0.0;
        public operator this(x:Double,y:Double):Double = x+y;
	}

    public static def main(args:Rail[String]) {
        new Boltzmann().boltzmann();
    }

    // TODO mutable Struct?
    // maps to g_fg
    public static final class LBDist {
        var actual:Double;
        var equilibrium:Double;
        var temp:Double;
    }

    // TODO mutable Struct?
    // maps to g_fld
    public static final class LBProperties {
        public var density:Double;
        /** momentum x */
        public var px:Double;
        /** momentum y */
        public var py:Double;
        public var pressure:Double;
        /** relaxation */
        public var t_rho:Double;

        public def this() { }
        public def this(density:Double, px:Double, py:Double, pressure:Double, t_rho:Double) {
            this.density = density;
            this.px = px;
            this.py = py;
            this.pressure = pressure;
            this.t_rho = t_rho;
        }
    }
}

