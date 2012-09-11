#include <stdlib.h>
#include <sys/time.h> 
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <papi.h>
#include "ElectrostaticDirect.h"
#include "Atom.h"

//#define USEPAPI 1

using namespace std;

static long microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}

static void printerror(const char *file, int line, const char *call, int retval);

static Atom* setup() {
	atoms = new Atom[atomsPerPlace+1];
    int gridSize = (int)(ceil(cbrt(numAtoms)));
    int gridPoint = myStart; // running total of assigned grid points
	for (int i=0; i<myNumAtoms; i++) {
		int gridX = gridPoint / (gridSize * gridSize);
        int gridY = (gridPoint - (gridX * gridSize * gridSize)) / gridSize;
        int gridZ = gridPoint - (gridX * gridSize * gridSize) - (gridY * gridSize);
        double x = (gridX + 0.5 + randomNoise()) * (SIZE / gridSize);
        double y = (gridY + 0.5 + randomNoise()) * (SIZE / gridSize);
        double z = (gridZ + 0.5 + randomNoise()) * (SIZE / gridSize);
        atoms[i].centre.i = x;
		atoms[i].centre.j = y;
		atoms[i].centre.k = z;
		atoms[i].charge = gridPoint%2==0?1.0:-1.0;
        gridPoint++;
	}
    return atoms;
}

static double randomNoise() {
	return 0.0;
}

static double getEnergy() {
#ifdef USEPAPI
    int EventSet=PAPI_NULL;
    int retval = PAPI_library_init(PAPI_VER_CURRENT);

    if (retval != PAPI_VER_CURRENT && retval > 0) {
        fprintf(stderr,"PAPI library version mismatch!\n");
        exit(1);
    }

    if (retval < 0) {
        fprintf(stderr, "Initialization error!\n");
        exit(1);
    }
    PAPI_create_eventset(&EventSet);
    PAPI_add_event(EventSet, PAPI_TOT_INS);
    PAPI_add_event(EventSet, PAPI_FP_INS);
    PAPI_add_event(EventSet, PAPI_FP_OPS);
    long long *totals = new long long[3];
    bzero(totals, 3*sizeof(long long));
    long long cycles = -PAPI_get_real_cyc();
    if ((retval=PAPI_start(EventSet)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_start", retval);
#endif

    double energy = 0.0;

    Atom* otherAtoms = new Atom[atomsPerPlace+1];
    Atom* nextAtoms = new Atom[atomsPerPlace+1];

    MPI_Request sendAtoms;
    int firstTarget = (myRank+1)%nTasks;
    MPI_Request recvAtoms;
    int firstSource = (myRank+nTasks-1)%nTasks;
    if (nTasks > 1) {
        // asynchronously send first set of atoms to next task
        MPI_Isend(atoms, (atomsPerPlace+1)*7, MPI_DOUBLE, firstTarget, 1, MPI_COMM_WORLD, &sendAtoms);

        // asynchronously receive first set of atoms from previous task
        MPI_Irecv(otherAtoms, (atomsPerPlace+1)*7, MPI_DOUBLE, firstSource, 1, MPI_COMM_WORLD, &recvAtoms);
    }

    // energy for all interactions within this place
    for (int i = 0; i < myNumAtoms; i++) {
        Atom atomI = atoms[i];
        double ix = atomI.centre.i;
        double iy = atomI.centre.j;
        double iz = atomI.centre.k;
        double qi = atomI.charge;

        double fi = 0.0;
        double fj = 0.0;
        double fk = 0.0;
        for (int j = 0; j < i; j++) {
            Atom atomJ = atoms[j];
            double dx = ix - atomJ.centre.i;
            double dy = iy - atomJ.centre.j;
            double dz = iz - atomJ.centre.k;

            double r2 = dx * dx + dy * dy + dz * dz;
            double invR2 = 1.0 / r2;
            double chargeInt = (qi * atomJ.charge);
            double invR = sqrt(invR2);

            double e = chargeInt * invR;
            double forceScale =  e * invR2;
            energy += 2.0 * e;

            double fx = forceScale * dx;
            double fy = forceScale * dy;
            double fz = forceScale * dz;

            fi += fx;
            fj += fy;
            fk += fz;

            atomJ.force.i -= fx;
            atomJ.force.j -= fy;
            atomJ.force.k -= fz;
        }
        atomI.force.i = fi;
        atomI.force.j = fj;
        atomI.force.k = fk;
    }

    MPI_Status ignore;
    if (nTasks > 1) {
        MPI_Wait(&sendAtoms, &ignore);
        MPI_Wait(&recvAtoms, &ignore);
    }

    int p = firstSource;

    for (int jump = 2; jump <= nTasks; jump++) {
        int source = (myRank+nTasks-jump)%nTasks;
        if (jump < nTasks) {
            int target = (myRank+jump)%nTasks;
            MPI_Isend(atoms, (atomsPerPlace+1)*7, MPI_DOUBLE, target, 2, MPI_COMM_WORLD, &sendAtoms);
            MPI_Irecv(nextAtoms, (atomsPerPlace+1)*7, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &recvAtoms);
        }

        // energy for all interactions with other atoms at other place
        int otherNumAtoms = p < leftOver ? atomsPerPlace+1 : atomsPerPlace;
        for (int i = 0; i < otherNumAtoms; i++) {
            for (int j = 0; j < myNumAtoms; j++) {
                double xDist = otherAtoms[i].centre.i - atoms[j].centre.i;
                double yDist = otherAtoms[i].centre.j - atoms[j].centre.j;
                double zDist = otherAtoms[i].centre.k - atoms[j].centre.k;

                double r2 = xDist * xDist + yDist * yDist + zDist * zDist;
                double invR2 = 1.0 / r2;
                double invR = sqrt(invR2);

                double e = (otherAtoms[i].charge * atoms[j].charge) * invR;
                energy += e;
                atoms[i].force.i += e * invR2 * xDist;
                atoms[i].force.j += e * invR2 * yDist;
                atoms[i].force.k += e * invR2 * zDist;
            }
        }

        if (jump < nTasks) {
            // finish receive of next set of atoms and swap working sets
            MPI_Wait(&recvAtoms, &ignore);
            MPI_Wait(&sendAtoms, &ignore);
            Atom* tmp = otherAtoms;
            otherAtoms = nextAtoms;
            nextAtoms = tmp;
            p = source;
        }
    }
/*        
    for (int p = 0; p < nTasks; p++) {
        // broadcast atoms from place p to other places
        if (p == myRank) {
            memcpy(otherAtoms, atoms, sizeof(Atom)*(atomsPerPlace+1));
        }
        MPI_Bcast(otherAtoms, (atomsPerPlace+1)*4, MPI_DOUBLE, p, MPI_COMM_WORLD);

        if (p != myRank) {
            // energy for all interactions with other atoms at other place
            int otherNumAtoms = p < leftOver ? astatic void printerror(const char *file, int line, const char *call, int retval);tomsPerPlace+1 : atomsPerPlace;
            for (int i = 0; i < otherNumAtoms; i++) {
                for (int j = 0; j < myNumAtoms; j++) {
                    double xDist = otherAtoms[i].centre.i - atoms[j].centre.i;
                    double yDist = otherAtoms[i].centre.j - atoms[j].centre.j;
                    double zDist = otherAtoms[i].centre.k - atoms[j].centre.k;

                    double distance = sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

                    energy += (otherAtoms[i].charge * atoms[j].charge) / distance;
                }
            }
        } else {
            // energy for all interactions within this place
            for (int i = 0; i < myNumAtoms; i++) {
                for (int j = 0; j < i; j++) {
                    double xDist = atoms[i].centre.i -atoms[j].centre.i;
                    double yDist = atoms[i].centre.j -atoms[j].centre.j;
                    double zDist = atoms[i].centre.k -atoms[j].centre.k;

                    double distance = sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

                    energy += 2.0 * (atoms[i].charge * atoms[j].charge) / distance;
                }
            }
        }
    }
*/

    energy /= 2.0;

    double totalEnergy;
    MPI_Reduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef USEPAPI
    cycles += PAPI_get_real_cyc();
    PAPI_stop(EventSet, totals);
    printf("cycles: %16lld total ins: %16lld FP ins: %16lld FLOPS: %16lld\n", cycles, totals[0], totals[1], totals[2]);
    printf("FLOPS/cycle %f\n", (double)totals[1] / cycles);
#endif
    
    if (myRank == 0) {
        printf("energy = %.5f\n", totalEnergy);
    }

    return energy;
}

static void printerror(const char *file, int line, const char *call, int retval) {
    printf("%s\tFAILED\nLine # %d\n", file, line);
    if ( retval == PAPI_ESYS ) {
        char buf[128];
        memset( buf, '\0', sizeof(buf) );
        sprintf(buf, "System error in %s:", call );
        perror(buf);
    }
    else if ( retval > 0 ) {
        printf("Error calculating: %s retval %d\n", call, retval );
    }
    else {
        char errstring[PAPI_MAX_STR_LEN];
        PAPI_perror(retval, errstring, PAPI_MAX_STR_LEN );
        printf("Error in %s: %s\n", call, errstring );
    }
    printf("\n");
    exit(1);
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

	numAtoms = 30000;
	if (argc > 1) {
		char* end;
		numAtoms = strtol(argv[1], &end, 10);
	}

    atomsPerPlace = numAtoms / nTasks;
    leftOver = numAtoms % nTasks;
    myNumAtoms = myRank < leftOver ? atomsPerPlace+1 : atomsPerPlace;
    myStart = myRank < leftOver ? (atomsPerPlace+1) * myRank : leftOver + atomsPerPlace*myRank;
    myEnd = myStart + myNumAtoms - 1;

    Atom* atoms = setup();
    long start = microTime();
	double energy = getEnergy();
    long stop = microTime();
    if (myRank == 0) {
        printf("%.3g s\n", (stop-start) / 1.0e6);
    }

    MPI_Finalize();
    return 0;
}
