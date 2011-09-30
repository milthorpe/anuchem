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
import x10.io.File;
import x10.io.FileWriter;
import x10.io.Printer;

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
    public static SIZE = new Array[Int](2, NSIZE);
    public static GHOST_WIDTH = 1;

    public static NSTEPS = 5000;

    public static VISCOSITY = 1.0;

    /** time increment*/
    public static DELTA_T = 1.0;

    /** the dimensions of the square cavity */
    public static XMAX = 256.0;

    /** the distance between lattice points */
    public static DELTA_X = XMAX/((NSIZE-1) as Double);

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
    val latticeRegion = 0..(NSIZE+1) * 0..(NSIZE+1);

    /** The distributed of the full lattice. */
    val latticeDist = Dist.makeBlockBlock(latticeRegion, 0, 1);

    /** The region of interior points of the lattice. */
    val interiorRegion:Region(2){rect};

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

        val cspd = Math.sqrt(2.0)*DELTA_X/DELTA_T;
        Console.OUT.println("cspd = " + cspd);
        this.cspd = cspd;
        val tau_rho = 6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5;
        Console.OUT.println("Value of tau_rho is " + tau_rho);

        Console.OUT.println("latticeDist = " + latticeDist);

        // initialize boundary conditions
        val boundaryCond = ([i,j]:Point(2)) => { (i==1 || i==(SIZE(0)) || j==1) ? 1 : (j==(SIZE(1))) ? 2 : 0 };
        val boundary = DistArray.make[Int](latticeDist, boundaryCond, GHOST_WIDTH);
        boundary.updateGhosts();
        this.boundary = boundary;

        val interiorRegion = 1..NSIZE * 1..NSIZE;

        // initialize distribution of density and velocities
        val lbProps = DistArray.make[LBProperties](latticeDist, 
            (p:Point) => new LBProperties(
                RHO0, // + 0.0*cos(2.0d00*pi*dble(j-1) / dble(size(2)-1))                                    
                ux0, uy0,
                RHO0 * RGAS * TMPRTR0 / (1.0 - B_VDW*RHO0) - A_VDW*RHO0*RHO0,
                6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5),
            GHOST_WIDTH);
        this.lbProps = lbProps;
        this.interiorRegion = interiorRegion;
        rtot = RHO0*lbProps.region.size();

        initializeLatticeParameters();

        val lbDist = Dist.makeBlockBlock(latticeRegion * 0..8, 0, 1);
        Console.OUT.println("lbDist = " + lbDist);
        val lb = DistArray.make[LBDist](lbDist, (Point) => new LBDist(), GHOST_WIDTH);
        this.lb = lb;
    }

    private def printLb() {
        for (place in lb.dist.places()) at (place) {
            val ffa = new Array[Double](9);
            val lbsLocal = lb.getLocalPortion();
            val lbHere = lb.dist(here) as Region(3){rect};
            for (ii in lbHere.min(0)..lbHere.max(0)) {
                for (jj in lbHere.min(1)..lbHere.max(1)) {
                    for (i in 0..8) {
                        Console.OUT.printf("%12.4f", lbsLocal(ii,jj,i).actual);
                    }
                    Console.OUT.println();
                }
                Console.OUT.println();
            }
        }
    }

    private def printLbProps() {
        for (place in lb.dist.places()) at (place) {
            val ffa = new Array[Double](9);
            val lbPropsLocal = lbProps.getLocalPortion();
            val latticeRegionHere = latticeDist(here) as Region(2){rect};
            for (i in latticeRegionHere.min(0)..latticeRegionHere.max(0)) {
                for (j in latticeRegionHere.min(1)..latticeRegionHere.max(1)) {
                    Console.OUT.println(lbPropsLocal(i,j).density);
                }
                Console.OUT.println();
            }
        }
    }
    

    public def boltzmann() {
        // initialize equilibrium and actual distributions
        equilibrium();
        finish ateach([ii,jj] in latticeDist) {
            for (i in 0..8) {
                lb(ii,jj,i).actual = lb(ii,jj,i).equilibrium;
            }
        }
        properties();
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
        hash(0,0)   = 0;
        hash(1,0)   = 1;
        hash(-1,0)  = 2;
        hash(0,1)   = 3;
        hash(0,-1)  = 4;
        hash(1,1)   = 5;
        hash(-1,-1) = 6;
        hash(1,-1)  = 7;
        hash(-1,1)  = 8;
        return hash;
    }

    private static def initIHash(): Array[Int](2){rect} {
        val ihash = new Array[Int](-1..1 * -1..1);
        ihash(0,0)   = 0;
        ihash(1,0)   = 2;
        ihash(-1,0)  = 1;
        ihash(0,1)   = 4;
        ihash(0,-1)  = 3;
        ihash(1,1)   = 6;
        ihash(-1,-1) = 5;
        ihash(1,-1)  = 8;
        ihash(-1,1)  = 7;
        return ihash;
    }

    private static def initEi(): Array[Int](2){rect,zeroBased} {
        val ei = new Array[Int](0..7 * 0..1);
        ei(0,0) =  1; ei(0,1) =  0;
        ei(1,0) = -1; ei(1,1) =  0;
        ei(2,0) =  0; ei(2,1) =  1;
        ei(3,0) =  0; ei(3,1) = -1;
        ei(4,0) =  1; ei(4,1) =  1;
        ei(5,0) = -1; ei(5,1) = -1;
        ei(6,0) =  1; ei(6,1) = -1;
        ei(7,0) = -1; ei(7,1) =  1;
        return ei;
    }

    /** Update values of density, momentum and pressure */
    private def properties() {
        val lb = this.lb; // TODO shouldn't be necessary XTENLANG-1913
        val latticeDist = this.latticeDist; // TODO shouldn't be necessary XTENLANG-1913
        val interiorRegion = this.interiorRegion; // TODO shouldn't be necessary XTENLANG-1913
        val lbProps = this.lbProps; // TODO shouldn't be necessary XTENLANG-1913

        val cspd2 = cspd / Math.sqrt(2.0);
		val rtot = new Accumulator[Double](SumReducer());
		val uxtot = new Accumulator[Double](SumReducer());
		val uytot = new Accumulator[Double](SumReducer());
        finish ateach (p in Dist.makeUnique()) {
            var myUxtot:Double = 0.0;
            var myUytot:Double = 0.0;
            var myRtot:Double = 0.0;
            val lbLocal = lb.getLocalPortion();
            val lbPropsLocal = lbProps.getLocalPortion();
            val interiorRegionHere = (latticeDist(here) && interiorRegion) as Region(2){rect};
            for([ii,jj] in interiorRegionHere) {
                val props = lbPropsLocal(ii,jj);
                props.density = 0.0;
                props.px = 0.0;
                props.py = 0.0;

                // evaluate density and momentum
                for (i in 0..8) {
                    val j = i-1;
                    val ex:Double;
                    val ey:Double;
                    if (j<0) {
                        ex = 0.0;
                        ey = 0.0;
                    } else {
                        ex = cspd2*ei(j,0);
                        ey = cspd2*ei(j,1);
                    }
                    val dd = lbLocal(ii,jj,i).actual;
                    props.density += dd;
                    props.px += ex * dd;
                    props.py += ey * dd;
                }
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

                myUxtot += props.px;
                myUytot += props.py;

                // evaluate pressure (eq. 4.9)
                val rho = props.density;
                props.pressure = rho*RGAS*TMPRTR0/(1.0-B_VDW*rho) - A_VDW*rho*rho;
                myRtot += rho;
                
            }
            uxtot <- myUxtot;
            uytot <- myUytot;
            rtot <- myRtot;
        }
        this.rtot = rtot();
        this.uxtot = uxtot();
        this.uytot = uytot();
    }

    /** Update equilibrium distributions */
    private def equilibrium() {
        val lb = this.lb; // TODO shouldn't be necessary XTENLANG-1913
        val latticeDist = this.latticeDist; // TODO shouldn't be necessary XTENLANG-1913
        val lbProps = this.lbProps; // TODO shouldn't be necessary XTENLANG-1913

        val rdim = 4.0;
        val b = 24.0;
        val c2 = cspd*cspd;
        val cspdi = Math.sqrt(2.0)/cspd;
        val cspd2 = 1.0/cspdi;

        finish ateach (p in Dist.makeUnique()) {
            val ffa = new Array[Double](9);
            val lbLocal = lb.getLocalPortion();
            val lbPropsLocal = lbProps.getLocalPortion();
            val latticeRegionHere = latticeDist(here) as Region(2){rect};
            for([ii,jj] in latticeRegionHere) {
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
                    val props = lbPropsLocal(ii,jj);
                    val aa = (rdim/(b*c2)) * props.pressure;
                    val aa0 = props.density - b*aa;
                    if (j == 0) {
                        ffa(j) = aa0 + 4.0*aa;
                    } else if (j < 5) {
                        ffa(j) = 4.0*aa;
                    } else {
                        ffa(j) = aa;
                    }

                    lbLocal(ii,jj,j).equilibrium = ffa(j) 
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

    /** Advance simulation one timestep */
    private def timestep() {
        timer.start(0);

        val lb = this.lb; // TODO shouldn't be necessary XTENLANG-1913
        val latticeDist = this.latticeDist; // TODO shouldn't be necessary XTENLANG-1913
        val boundary = this.boundary; // TODO shouldn't be necessary XTENLANG-1913
        val lbProps = this.lbProps; // TODO shouldn't be necessary XTENLANG-1913

        // make backup copies of distribution and update ghost cells
        //finish ateach ([ii,jj] in latticeDist | interiorRegion) {
        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbRegionHere = lb.dist(here) as Region(3){rect};
            for ([ii,jj,i] in lbRegionHere) {
                lbLocal(ii,jj,i).temp = lbLocal(ii,jj,i).actual;
            }
        }

        timer.start(1);
        lb.updateGhosts();
        timer.stop(1);

        // Perform streaming operation
        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion() as Array[LBDist](3){rect};
            val interiorRegionHere = (latticeDist(here) && interiorRegion) as Region(2){rect};
            val patch = new Array[Double](-1..1 * -1..1 * 0..8);
            for ([ii,jj] in interiorRegionHere) {
                if (boundary(ii,jj) == 0) {
                    for (i in 1..8) {
                        val k = i-1;
                        val ix = ei(k,0); //Math.round(ei(k,0)) as Int;
                        val iy = ei(k,1); //Math.round(ei(k,1)) as Int;
                        lbLocal(ii,jj,i).actual = lbLocal(ii-ix,jj-iy,i).temp;
                    }
                } else {
                    getPatch(ii, jj, patch, lbLocal); // maps to fgp
                    for (i in 1..8) {
                        val k = i-1;
                        val ix = ei(k,0); //Math.round(ei(k,0)) as Int;
                        val iy = ei(k,1); //Math.round(ei(k,1)) as Int;
                        lbLocal(ii,jj,i).actual = patch(-ix,-iy,i);
                    }   
                }
            }
        }

        properties();

        // perform relaxation
        equilibrium();

        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbPropsLocal = lbProps.getLocalPortion();
            val interiorRegionHere = (lb.dist(here) && (interiorRegion * 0..8)) as Region(3){rect};
            for ([ii,jj,i] in interiorRegionHere) {
                val t_rho = lbPropsLocal(ii,jj).t_rho;
                if (t_rho > 0.0) {
                    val t = lbLocal(ii,jj,i).temp;
                    lbLocal(ii,jj,i).actual = t - (t - lbLocal(ii,jj,i).equilibrium) / t_rho;
                }
                // TODO boundary constant condition is commented out in GA code
    //            if (bc(ii,jj) == 2) {
    //                lbLocal(ii,jj,i).actual = lbLocal(ii,jj,i).equilibrium;
    //            }
            }
        }
        timer.stop(0);
    }

    /** Handle cells at boundary */
    private def getPatch(ii:Int, jj:Int, patch:Array[Double](3){rect}, lbLocal:Array[LBDist](3){rect}) {
        // Check values of neighboring cells
        val mask = new Array[Int](-1..1 * -1..1, ([k,l]:Point(2)) =>
                (latticeRegion.contains(ii+k, jj+l) && boundary(ii+k, jj+l) == 0) ? 0 : 2);

        // Determine if cells in mask represent interior cells
        // or are on the boundary
        for ([k,l] in mask) {
            if (mask(k,l) != 0 && (k!=0 || l!=0)) {
                for (kk in Math.max(k-1,-1)..Math.min(k+1,1)) {
                    for (ll in Math.max(l-1,-1)..Math.min(l+1,1)) {
                        if (mask(kk,ll) == 0) mask(k,l) = 1;
                    }
                }
            }
        }

            
        // Evaluate distribution in boundary patch
        val bval = boundary(ii,jj);

        if (bval == 1) {
            // Apply simple bounce back condition
            for ([k,l] in mask) {
                if (mask(k,l) == 2) {
                    patch(k,l,ihash(k,l)) = lbLocal(ii,jj,hash(k,l)).temp;
                } else {
                    patch(k,l,ihash(k,l)) = lbLocal(ii+k,jj+l,ihash(k,l)).temp;
                }
            }

        } else if (bval == 2) {
            // Apply constant velocity boundary condition to flat upper boundary
            val cspd2 = cspd / Math.sqrt(2.0);
            var rho:Double = 0.0;
            var ux:Double = 0.0;
            var uy:Double = 0.0;
            for ([k,l] in mask) {
                if (mask(k,l) == 2) {
                    patch(k,l,ihash(k,l)) = lbLocal(ii,jj,hash(k,l)).temp;
                } else {
                    patch(k,l,ihash(k,l)) = lbLocal(ii+k,jj+l,ihash(k,l)).temp;
                }
                rho = rho + patch(k,l,ihash(k,l));
                ux = ux - cspd2 * k * patch(k,l,ihash(k,l));
                uy = uy - cspd2 * l * patch(k,l,ihash(k,l));
            }

            ux = ux / rho;
            uy = uy / rho;
            ux = UXBC - ux;

            // Add corrections needed to adjust for velocity mismatch
            patch(1,1,ihash(1,1)) = patch(1,1,ihash(1,1)) - 0.5*rho*ux/cspd2;
            patch(-1,1,ihash(-1,1)) = patch(-1,1,ihash(-1,1)) + 0.5*rho*ux/cspd2;
        
        } else if (bval == 3) {
            // Apply constant velocity boundary condition
            val cspd2 = Math.sqrt(2.0)*cspd;
            var rho:Double = 0.0;
            var ux:Double = 0.0;
            var uy:Double = 0.0;
            for ([k,l] in mask) {
                if (mask(k,l) == 2 && (k!=0 || l!=0)) {
                  patch(k,l,ihash(k,l)) = lbLocal(ii,jj,hash(k,l)).actual;
                  rho = rho + lbLocal(ii,jj,hash(k,l)).actual;
                  ux = ux + cspd2*k*lbLocal(ii,jj,hash(k,l)).actual;
                  uy = uy + cspd2*l*lbLocal(ii,jj,hash(k,l)).actual;
                } else {
                  patch(k,l,ihash(k,l)) = lbLocal(ii+k,jj+l,ihash(k,l)).actual; // TODO check: GA has ihash(i,l)
                  rho = rho + lbLocal(ii+k,jj+l,ihash(k,l)).actual;
                  ux = ux + cspd2*k*lbLocal(ii+k,jj+l,ihash(k,l)).actual;
                  uy = uy + cspd2*l*lbLocal(ii+k,jj+l,ihash(k,l)).actual;
                }
            }

            // Determine value of correction needed to get specified final velocity
            // TODO these are never used in GA code
            //val dux = bcpar(2,1) - ux;
            //val duy = bcpar(2,2) - uy;
        }
    }

    /** Evaluate the vorticity of velocity field */
    private def vorticity() {
        lbProps.updateGhosts();
        finish ateach (p in Dist.makeUnique()) {
            val lbLocal = lb.getLocalPortion();
            val lbPropsLocal = lbProps.getLocalPortion();
            val lbInterior = (latticeDist.region && interiorRegion) as Region(2){rect};
            for ([ii,jj] in lbInterior) {
                val propsThis = lbPropsLocal(ii,jj);
                var drni:Double = 0.0;
                var drxi:Double = 0.0;
                var dryi:Double = 0.0;
                var drho_x:Double = 0.0;
                var drho_y:Double = 0.0;
                var dux_x:Double = 0.0;
                var dux_y:Double = 0.0;
                var duy_x:Double = 0.0;
                var duy_y:Double = 0.0;
                for (j in -1..1) {
                    val dy = j * DELTA_X;
                    val dyi:Double;
                    if (dy == 0.0) {
                        dyi = 0.0;
                    } else {
                        dyi = 1.0 / dy;
                    }
                    for (i in -1..1) {
                        if (boundary(ii+i,jj+j) == 0 && (i!=0||j!=0)) {
                            val dx = i * DELTA_X;
                            val dxi:Double;
                            if (dx == 0.0) {
                                dxi = 0.0;
                            } else {
                                dxi = 1.0 / dx;
                            }
                            val dr = Math.sqrt(dx*dx + dy*dy);
                            val dri = 1.0/dr;
                            val adxi = Math.abs(dx*dri*dri);
                            val adyi = Math.abs(dy*dri*dri);
                            val propsOther = lbPropsLocal(ii+i,jj+j);
                            val rhot = propsOther.density;
                            val uxt = propsOther.px;
                            val uyt = propsOther.py;
                            drho_x += (rhot-propsThis.density)*dxi*adxi;
                            drho_y += (rhot-propsThis.density)*dyi*adyi;
                            dux_x += (uxt-propsThis.px)*dxi*adxi;
                            dux_y += (uxt-propsThis.px)*dyi*adyi;
                            duy_x += (uyt-propsThis.py)*dxi*adxi;
                            duy_y += (uyt-propsThis.py)*dyi*adyi;
                            drxi += adxi;
                            dryi += adyi;
                        }
                    }
                }

                drho_x /= drxi;
                drho_y /= dryi;
                dux_x /= drxi;
                dux_y /= dryi;
                duy_x /= drxi;
                duy_y /= dryi;
                propsThis.vorticity = duy_x - dux_y;
            }
        }
    }

    /** Print current value of fields to a file */
    private def printData() {
        val imax = NSIZE/10;
        val glo = new Array[Int](3, 1);
        val ghi = new Array[Int](3);
        ghi(0) = NSIZE;
        ghi(2) = 3;
        val bld = new Array[Int](2);
        bld(0) = NSIZE;
        bld(1) = 10;

        val MAXELEM=256;

        // Check dimensions to see if size needs to be reduced
        val inc1:Int;
        val inc2:Int;
        val icnt1:Int;
        val icnt2:Int;
        // TODO non-square array (use size(i))
        if (NSIZE > MAXELEM) {
            inc1 = SIZE(0)/MAXELEM;
            inc2 = SIZE(1)/MAXELEM;
            icnt1 = SIZE(0)/inc1;
            icnt2 = SIZE(1)/inc2;
        } else {
            inc2 = inc1 = 1;
            icnt1 = SIZE(0);
            icnt2 = SIZE(1);
        }


        val fil = new Printer(new FileWriter(new File("bltz.gmv")));
        fil.println("gmvinput ascii");
        fil.printf("nodes       -1 %10i%10i%10i\n", icnt1, icnt2, 1);

        val dx = inc1*XMAX/(SIZE(0)-1);
        val dy = inc1*XMAX/(SIZE(1)-1);
        for (i in 0..(icnt1-1)) {
            fil.printf("    %12.4f", i*dx);
            if (i%5 == 4) fil.println();
        }
        fil.println();
        for (i in 0..(icnt2-1)) {
            fil.printf("    %12.4f", i*dy);
            if (i%5 == 4) fil.println();
        }
        fil.println();
        fil.printf("    %12.4f", 0.0);
        fil.println("cells     0");
        fil.println("variable");
        fil.println("rho                           1");
        var jcnt:Int = 0;
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
        public var vorticity:Double;

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

