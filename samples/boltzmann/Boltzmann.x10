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
import x10.compiler.StackAllocate;
import x10.io.File;
import x10.io.FileWriter;
import x10.io.Printer;
import x10.regionarray.Array;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;
import x10.util.Team;

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
public class Boltzmann(nsize:Long, nsteps:Int) {
    public static GHOST_WIDTH = 1;

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

    public val size:Rail[Long];

    /** the distance between lattice points */
    public val delta_x:Double;

    /** The distribution of the full lattice. */
    val latticeDist:Dist(2);

    /** The lattice Boltzmann distribution functions at each site. maps to g_fg(..,..,1..9) */
    val current:DistArray[Double](3);

    /** The equilibrium distribution functions at each site. maps to g_fg(..,..,10..18) */
    val equilibrium:DistArray[Double](3){dist==current.dist};

    /** The previous values of distribution functions. maps to g_fg(..,..,19..27) */
    val previous:DistArray[Double](3){dist==current.dist};

    /** The properties (density, momentum etc.) at each site. maps to g_fld*/
    val lbProps:DistArray[LBProperties](2);

    /** The boundary condition data. maps to g_bc*/
    val boundary:DistArray[Long](2);

    /** Speed - magnitude of velocity for e_i */
    val cspd:Double;

    public val timer = PlaceLocalHandle.make[Timer](Place.places(), () => new Timer(2)); // maps to tstats

    val ffb = new Rail[Double](9);
    val ffc = new Rail[Double](9);
    val ffd = new Rail[Double](9);

    /** The velocity vectors e_i */
    val ei = initEi();

    static hash = initHash();
    static ihash = initIHash();

    public def this(nsize:Long, nsteps:Int) {
        property(nsize, nsteps);
        val size = new Rail[Long](2, nsize);
        this.size = size;
        val delta_x = XMAX/((nsize-1) as Double);
        this.delta_x = delta_x;

        val latticeRegion = Region.makeRectangular(0..(size(0)-1), 0..(size(1)-1));
        val latticeDist = Dist.makeBlockBlock(latticeRegion, 0, 1);
        this.latticeDist = latticeDist;

        val cspd = Math.sqrt(2.0)*delta_x/DELTA_T;
        Console.OUT.println("cspd = " + cspd);
        this.cspd = cspd;
        val tau_rho = 6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5;
        Console.OUT.println("Value of tau_rho is " + tau_rho);

        Console.OUT.println("latticeDist = " + latticeDist);

        // initialize boundary conditions
        val boundaryCond = ([i,j]:Point(2)) => { (i==0 || i==(size(0)-1) || j==0) ? 1 : (j==(size(1)-1)) ? 2 : 0 };
        val boundary = DistArray.make[Long](latticeDist, boundaryCond, GHOST_WIDTH, false);
        boundary.updateGhosts();
        this.boundary = boundary;

        // initialize distribution of density and velocities
        val ux0 = 0.0;
        val uy0 = 0.0;
        val lbProps = DistArray.make[LBProperties](latticeDist, 
            (p:Point) => new LBProperties(
                RHO0, // + 0.0*cos(2.0d00*pi*dble(j-1) / dble(size(2)-1))                                    
                ux0, uy0,
                RHO0 * RGAS * TMPRTR0 / (1.0 - B_VDW*RHO0) - A_VDW*RHO0*RHO0,
                6.0 * VISCOSITY / (cspd*cspd*DELTA_T) + 0.5),
            GHOST_WIDTH,
            false);
        this.lbProps = lbProps;
        val rtot = RHO0*lbProps.region.size();

        initializeLatticeParameters();

        val lbDist = Dist.makeBlockBlock(latticeRegion * Region.makeRectangular(0..8), 0, 1);
        Console.OUT.println("lbDist = " + lbDist);
        this.current = DistArray.make[Double](lbDist, 0, false);
        this.equilibrium = DistArray.make[Double](lbDist, 0, false);
        this.previous = DistArray.make[Double](lbDist, GHOST_WIDTH, false);
    }

    private def printDistributions() {
        Console.OUT.println("new timestep");
        for (place in current.dist.places()) at (place) {
            val currentLocal = current.getLocalPortion();
            val lbHere = current.dist(here) as Region(3){rect};
            for (jj in lbHere.min(1)..lbHere.max(1)) {
                for (ii in lbHere.min(0)..lbHere.max(0)) {
                    Console.OUT.print("    ");
                    for (i in 0..8) {
                        Console.OUT.printf("%12.4f", currentLocal(ii,jj,i));
                    }
                    Console.OUT.println();
                }
                Console.OUT.println();
            }
        }
    }

    private def printLbProps() {
        val lbProps = this.lbProps;
        val latticeDist = this.latticeDist;
        for (place in lbProps.dist.places()) at (place) {
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
        // initialize distributions
		val stats = new Accumulator[LBStatistics](LBStatistics.SumReducer());

        finish ateach (p in Dist.makeUnique()) {
            val currentLocal = current.getLocalPortion() as Array[Double](3){rect};
            val equilibriumLocal = equilibrium.getLocalPortion() as Array[Double](3){rect};
            val previousLocal = previous.getLocalPortion() as Array[Double](3){rect};
            val lbPropsLocal = lbProps.getLocalPortion() as Array[LBProperties](2){rect};
            val latticeRegionLocal = latticeDist(here) as Region(2){rect};

            updateEquilibrium(latticeRegionLocal, equilibriumLocal, lbPropsLocal);
            for ([ii,jj] in latticeRegionLocal) {
                for (i in 0..8) {
                    previousLocal(ii,jj,i) = currentLocal(ii,jj,i) = equilibriumLocal(ii,jj,i);
                }
            }
            stats <- updateProperties(latticeRegionLocal, currentLocal, lbPropsLocal);
        }
        val initialStats = stats();
        Console.OUT.println("Total mass is " + initialStats.rtot);
        Console.OUT.println("Total x-momentum is " + initialStats.uxtot);
        Console.OUT.println("Total y-momentum is " + initialStats.uytot);

        finish ateach (p in Dist.makeUnique()) {
            for (i in 1..nsteps) {
                val myStats = timestep();
                if (i % 200 == 0) {
                    val currentStats = myStats.teamReduce();
                    if (here == Place.FIRST_PLACE) {
                        Console.OUT.println("Completed step " + i);
                        Console.OUT.println("Total density is " + currentStats.rtot);
                        Console.OUT.println("Total x-momentum is " + currentStats.uxtot);
                        Console.OUT.println("Total y-momentum is " + currentStats.uytot);
                    }
                    vorticity();
                    printData();
                }
            }
            vorticity();
            printData();
        }
        Console.OUT.printf("Time in lattice Boltzmann updates :%12.4f\n", (timer().total(0) as Double) / 1e9);
        Console.OUT.printf("Time in ghost cell updates        :%12.4f\n", (timer().total(1) as Double) / 1e9);
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

    private static def initHash():Array[Int](2){rect} {
        val hash = new Array[Int](MASK_REGION);
        hash(0,0)   = 0n;
        hash(1,0)   = 1n;
        hash(-1,0)  = 2n;
        hash(0,1)   = 3n;
        hash(0,-1)  = 4n;
        hash(1,1)   = 5n;
        hash(-1,-1) = 6n;
        hash(1,-1)  = 7n;
        hash(-1,1)  = 8n;
        return hash;
    }

    private static def initIHash():Array[Int](2){rect} {
        val ihash = new Array[Int](MASK_REGION);
        ihash(0,0)   = 0n;
        ihash(1,0)   = 2n;
        ihash(-1,0)  = 1n;
        ihash(0,1)   = 4n;
        ihash(0,-1)  = 3n;
        ihash(1,1)   = 6n;
        ihash(-1,-1) = 5n;
        ihash(1,-1)  = 8n;
        ihash(-1,1)  = 7n;
        return ihash;
    }

    static struct DisplacementVector(x:Long, y:Long) { }

    private static def initEi():Rail[DisplacementVector] {
        val ei = new Rail[DisplacementVector](9);
        ei(0) = DisplacementVector( 0, 0);
        ei(1) = DisplacementVector( 1, 0);
        ei(2) = DisplacementVector(-1, 0);
        ei(3) = DisplacementVector( 0, 1);
        ei(4) = DisplacementVector( 0,-1);
        ei(5) = DisplacementVector( 1, 1);
        ei(6) = DisplacementVector(-1,-1);
        ei(7) = DisplacementVector( 1,-1);
        ei(8) = DisplacementVector(-1, 1);
        return ei;
    }

    /** Update values of density, momentum and pressure for a local region */
    private def updateProperties(latticeRegionLocal:Region(2){rect},
                           currentLocal:Array[Double](3){rect}, 
                           lbPropsLocal:Array[LBProperties](2){rect}):LBStatistics {
        val cspd2 = cspd / Math.sqrt(2.0);
        var myUxtot:Double = 0.0;
        var myUytot:Double = 0.0;
        var myRtot:Double = 0.0;
        for([ii,jj] in latticeRegionLocal) {
            val props = lbPropsLocal(ii,jj);
            props.density = 0.0;
            props.px = 0.0;
            props.py = 0.0;

            // evaluate density and momentum
            for (i in 0..8) {
                val displacement = ei(i);
                val ex = cspd2*displacement.x;
                val ey = cspd2*displacement.y;
                val dd = currentLocal(ii,jj,i);
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

            myUxtot += props.px * props.density;
            myUytot += props.py * props.density;

            // evaluate pressure (eq. 4.9)
            val rho = props.density;
            props.pressure = rho*RGAS*TMPRTR0/(1.0-B_VDW*rho) - A_VDW*rho*rho;
            myRtot += rho;
            
        }
        return LBStatistics(myRtot, myUxtot, myUytot);
    }

    private def updateEquilibrium(latticeRegionLocal:Region(2){rect},
                           equilibriumLocal:Array[Double](3){rect}, 
                           lbPropsLocal:Array[LBProperties](2){rect}) {
        val rdim = 4.0;
        val b = 24.0;
        val c2 = cspd*cspd;
        val cspd2 = cspd/Math.sqrt(2.0);

        for([ii,jj] in latticeRegionLocal) {
            val props = lbPropsLocal(ii,jj);
            val aa = (rdim/(b*c2)) * props.pressure;
            val aa0 = props.density - b*aa;
            val px = props.px;
            val py = props.py;
            for (j in 0..8) {
                val ex = cspd2*ei(j).x;
                val ey = cspd2*ei(j).y;
                val ffa:Double;
                if (j == 0) {
                    ffa = aa0 + 4.0*aa;
                } else if (j < 5) {
                    ffa = 4.0*aa;
                } else {
                    ffa = aa;
                }

                equilibriumLocal(ii,jj,j) = ffa 
                   + props.density * ffb(j)*(px*ex + py*ey)
                   + ffc(j) * (px*px * ex*ex + 2.0*px*py*ex*ey + py*py * ey*ey)
                   + ffd(j) * (px*px + py*py);
            }
        }
    }

    /** Advance simulation one timestep; called at each place */
    private def timestep():LBStatistics {
        val localTimer = timer();
        localTimer.start(0);

        localTimer.start(1);
        previous.sendGhostsLocal();
        previous.waitForGhostsLocal();
        localTimer.stop(1);

        val currentLocal = current.getLocalPortion() as Array[Double](3){rect};
        val equilibriumLocal = equilibrium.getLocalPortion() as Array[Double](3){rect};
        val previousLocal = previous.getLocalPortion() as Array[Double](3){rect};
        val lbPropsLocal = lbProps.getLocalPortion() as Array[LBProperties](2){rect};
        val latticeRegionLocal = latticeDist(here) as Region(2){rect};
        val latticeRegion = latticeDist.region as Region(2){rect,zeroBased};

        // Perform streaming operation
        val boundaryLocal = boundary.getLocalPortion() as Array[Long](2){rect};
        val patch = new Array[Double](MASK_REGION * Region.makeRectangular(0..8)); // maps to fgp
        val mask = new Array[Long](MASK_REGION);
        for ([ii,jj] in latticeRegionLocal) {
            if (boundaryLocal(ii,jj) == 0) {
                for (i in 1..8) {
                    val ix = ei(i).x; //Math.round(ei(i).x) as Long;
                    val iy = ei(i).y; //Math.round(ei(i).y) as Long;
                    currentLocal(ii,jj,i) = previousLocal(ii-ix,jj-iy,i);
                }
            } else {
                getPatch(ii, jj, patch, mask, currentLocal, previousLocal, boundaryLocal, latticeRegion);
                for (i in 1..8) {
                    val ix = ei(i).x; //Math.round(ei(i).x) as Long;
                    val iy = ei(i).y; //Math.round(ei(i).y) as Long;
                    currentLocal(ii,jj,i) = patch(-ix,-iy,i);
                }   
            }
        }

        val localStats = updateProperties(latticeRegionLocal, currentLocal, lbPropsLocal);

        // perform relaxation
        updateEquilibrium(latticeRegionLocal, equilibriumLocal, lbPropsLocal);

        for ([ii,jj] in latticeRegionLocal) {
            val tRho = lbPropsLocal(ii,jj).t_rho;
            if (tRho > 0.0) {
                val invTRho = 1.0 / tRho;
                for (i in 0..8) {
                    val t = previousLocal(ii,jj,i);
                    currentLocal(ii,jj,i) = t - (t - equilibriumLocal(ii,jj,i)) * invTRho;
                    // make backup copy of distribution
                    previousLocal(ii,jj,i) = currentLocal(ii,jj,i);
                }
            }
            // TODO boundary constant condition is commented out in GA code
//            if (bc(ii,jj) == 2) {
//                currentLocal(ii,jj,i) = equilibriumLocal(ii,jj,i);
//            }
        }
        localTimer.stop(0);

        return localStats;
    }

    static MASK_REGION = Region.makeRectangular(-1..1, -1..1);

    /** Handle cells at boundary */
    private def getPatch(ii:Long, jj:Long, patch:Array[Double](3){rect}, mask:Array[Long](2){rect}, currentLocal:Array[Double](3){rect}, previousLocal:Array[Double](3){rect}, boundaryLocal:Array[Long](2){rect}, latticeRegion:Region(2){rect,zeroBased}) {
        // Check values of neighboring cells
        for (k in -1..1) {
            for (l in -1..1) {
                if (latticeRegion.contains(ii+k, jj+l) && boundaryLocal(ii+k, jj+l) == 0) {
                    mask(k,l) = 0;
                } else {
                    mask(k,l) = 2;
                }
            }
        }

        // Determine if cells in mask represent interior cells
        // or are on the boundary
        for (k in -1..1) {
            for (l in -1..1) {
                if ((k!=0 || l!=0) && mask(k,l) != 0) {
                    for (kk in Math.max(k-1,-1)..Math.min(k+1,1)) {
                        for (ll in Math.max(l-1,-1)..Math.min(l+1,1)) {
                            if (mask(kk,ll) == 0) mask(k,l) = 1;
                        }
                    }
                }
            }
        }
            
        // Evaluate distribution in boundary patch
        val bval = boundaryLocal(ii,jj);
        if (bval == 1) {
            // Apply simple bounce back condition
            for (k in -1..1) {
                for (l in -1..1) {
                    if (mask(k,l) == 2) {
                        patch(k,l,ihash(k,l)) = previousLocal(ii,jj,hash(k,l));
                    } else {
                        patch(k,l,ihash(k,l)) = previousLocal(ii+k,jj+l,ihash(k,l));
                    }
                }
            }

        } else if (bval == 2) {
            // Apply constant velocity boundary condition to flat upper boundary
            val cspd2 = cspd / Math.sqrt(2.0);
            var rho:Double = 0.0;
            var ux:Double = 0.0;
            var uy:Double = 0.0;
            for (k in -1..1) {
                for (l in -1..1) {
                    val p_kl:Double;
                    if (mask(k,l) == 2) {
                        p_kl = previousLocal(ii,jj,hash(k,l));
                    } else {
                        p_kl = previousLocal(ii+k,jj+l,ihash(k,l));
                    }
                    patch(k,l,ihash(k,l)) = p_kl;
                    rho += p_kl;
                    ux -= cspd2 * k * p_kl;
                    uy -= cspd2 * l * p_kl;
                }
            }

            ux /= rho;
            uy /= rho;
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
            for (k in -1..1) {
                for (l in -1..1) {
                    if (mask(k,l) == 2 && (k!=0 || l!=0)) {
                        val p_kl = currentLocal(ii,jj,hash(k,l));
                        patch(k,l,ihash(k,l)) = p_kl;
                        rho += p_kl;
                        ux += cspd2*k*p_kl;
                        uy += cspd2*l*p_kl;
                    } else {
                        val p_kl = currentLocal(ii+k,jj+l,ihash(k,l));
                        patch(k,l,ihash(k,l)) = p_kl; // TODO check: GA has ihash(i,l)
                        rho += p_kl;
                        ux += cspd2*k*p_kl;
                        uy += cspd2*l*p_kl;
                    }
                }
            }

            // Determine value of correction needed to get specified final velocity
            // TODO these are never used in GA code
            //val dux = bcpar(2,1) - ux;
            //val duy = bcpar(2,2) - uy;
        }
    }

    /** Evaluate the vorticity of velocity field; called at each place */
    private def vorticity() {
        val latticeDist = this.latticeDist; // TODO shouldn't be necessary XTENLANG-1913
        val lbProps = this.lbProps; // TODO shouldn't be necessary XTENLANG-1913

        lbProps.sendGhostsLocal();
        lbProps.waitForGhostsLocal();

        val lbPropsLocal = lbProps.getLocalPortion();
        val interiorRegion = Region.makeRectangular(1..(size(0)-2), 1..(size(1)-2));
        val interiorRegionHere = (interiorRegion && latticeDist(here)) as Region(2){rect};
        for ([ii,jj] in interiorRegionHere) {
            val props = lbPropsLocal(ii,jj);
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
                val dy = j * delta_x;
                val dyi:Double;
                if (dy == 0.0) {
                    dyi = 0.0;
                } else {
                    dyi = 1.0 / dy;
                }
                for (i in -1..1) {
                    if (boundary(ii+i,jj+j) == 0 && (i!=0||j!=0)) {
                        val dx = i * delta_x;
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
                        drho_x += (rhot-props.density)*dxi*adxi;
                        drho_y += (rhot-props.density)*dyi*adyi;
                        dux_x += (uxt-props.px)*dxi*adxi;
                        dux_y += (uxt-props.px)*dyi*adyi;
                        duy_x += (uyt-props.py)*dxi*adxi;
                        duy_y += (uyt-props.py)*dyi*adyi;
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
            props.vorticity = duy_x - dux_y;
        }
    }

    /** Print current value of fields to a file; called at each place */
    private def printData() {
        if (here == Place.FIRST_PLACE) {
            val imax = nsize/10;
            val glo = new Rail[Long](3, 1);
            val ghi = new Rail[Long](3);
            ghi(0) = nsize;
            ghi(2) = 3;
            val bld = new Rail[Long](2);
            bld(0) = nsize;
            bld(1) = 10;

            val MAXELEM=256;

            // Check dimensions to see if size needs to be reduced
            val inc1:Long;
            val inc2:Long;
            val icnt1:Long;
            val icnt2:Long;
            // TODO non-square array (use size(i))
            if (nsize > MAXELEM) {
                inc1 = size(0)/MAXELEM;
                inc2 = size(1)/MAXELEM;
                icnt1 = size(0)/inc1;
                icnt2 = size(1)/inc2;
            } else {
                inc2 = inc1 = 1;
                icnt1 = size(0);
                icnt2 = size(1);
            }

            val fil = new Printer(new FileWriter(new File("bltz.gmv")));
            fil.println("gmvinput ascii");
            fil.printf("nodes       -1%10i%10i%10i", icnt1, icnt2, 1);

            val dx = inc1*XMAX/(size(0)-1);
            val dy = inc1*XMAX/(size(1)-1);
            for (i in 0..(icnt1-1)) {
                if (i%5 == 0) fil.print("\n    ");
                fil.printf("%12.4f", i*dx);
            }
            for (i in 0..(icnt2-1)) {
                if (i%5 == 0) fil.print("\n    ");
                fil.printf("%12.4f", i*dy);
            }
            fil.printf("\n    %12.4f\n", 0.0);
            fil.println("cells     0");
            fil.println("variable");
            fil.println("rho                           1");

    /*
            var jcnt:Long = 0;
            for (i in 0..(imax-1)) {
                glo(2) = (i-1)*10 + 1;
                ghi(2) = i*10;
                if (ghi(2) > size(2)) ghi(2) = size(2);
                if (ghi(2) >= glo(2)) {
                    //call nga_get(g_fld, glo, ghi, buffer, bld)
                    for (k in 0..(ghi(2)-glo(2))) {
                        jcnt++;                      
                        if (((jcnt-1) % inc2) == 0) {
                            for (j in 0..(size(1)-1) {
                                if (i%5 == 0) fil.print("\n    ");
                                fil.printf("%12.4f", lbProps(i,j));
                            }
                        }
                    }
                }
            }
    */
            fil.println("ux                            1");
        }
    }

    public static def main(args:Rail[String]) {
        var nsize:Long = 129;
        var nsteps:Int = 5000n;
        if (args.size > 0) {
            nsize = Long.parseLong(args(0));
            if (args.size > 1) {
                nsteps = Int.parseInt(args(1));
            }
        }
        new Boltzmann(nsize, nsteps).boltzmann();
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

    /** 
     * This struct holds aggregate properties for a local portion held at a
     * place, or for the system as a whole.
     */
    public static struct LBStatistics(
        /* total density */
        rtot:Double,
        /* total x momentum */
        uxtot:Double,
        /* total y momentum */
        uytot:Double
    ) {
        public def this(rtot:Double, uxtot:Double, uytot:Double) {
            property(rtot, uxtot, uytot);
        }
    
        public static struct SumReducer implements Reducible[LBStatistics] {
            public def zero():LBStatistics = LBStatistics(0.0, 0.0, 0.0);
            public operator this(a:LBStatistics,b:LBStatistics):LBStatistics {
                return LBStatistics(a.rtot  + b.rtot,
                                    a.uxtot + b.uxtot,
                                    a.uytot + b.uytot);
            }
	    }

        public def teamReduce():LBStatistics {
            val stats = new Rail[Double](3);
            stats(0) = rtot;
            stats(1) = uxtot;
            stats(2) = uytot;
            Team.WORLD.allreduce[Double](stats, 0, stats, 0, 3, Team.ADD);
            return new LBStatistics(stats(0), stats(1), stats(2));
        }
    }
}

