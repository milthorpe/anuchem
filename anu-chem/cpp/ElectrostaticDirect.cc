#include <stdlib.h>
#include <sys/time.h> 
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "ElectrostaticDirect.h"
#include "Atom.h"

using namespace std;

static long microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}


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
    double energy = 0.0;

    Atom* otherAtoms = new Atom[atomsPerPlace+1];
    Atom* nextAtoms = new Atom[atomsPerPlace+1];

    MPI_Request sendAtoms;
    int firstTarget = (myRank+1)%nTasks;
    MPI_Request recvAtoms;
    int firstSource = (myRank+nTasks-1)%nTasks;
    if (nTasks > 1) {
        // asynchronously send first set of atoms to next task
        MPI_Isend(atoms, (atomsPerPlace+1)*4, MPI_DOUBLE, firstTarget, 1, MPI_COMM_WORLD, &sendAtoms);

        // asynchronously receive first set of atoms from previous task
        MPI_Irecv(otherAtoms, (atomsPerPlace+1)*4, MPI_DOUBLE, firstSource, 1, MPI_COMM_WORLD, &recvAtoms);
    }

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
            MPI_Isend(atoms, (atomsPerPlace+1)*4, MPI_DOUBLE, target, 2, MPI_COMM_WORLD, &sendAtoms);
            MPI_Irecv(nextAtoms, (atomsPerPlace+1)*4, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &recvAtoms);
        }

        // energy for all interactions with other atoms at other place
        int otherNumAtoms = p < leftOver ? atomsPerPlace+1 : atomsPerPlace;
        for (int i = 0; i < otherNumAtoms; i++) {
            for (int j = 0; j < myNumAtoms; j++) {
                double xDist = otherAtoms[i].centre.i - atoms[j].centre.i;
                double yDist = otherAtoms[i].centre.j - atoms[j].centre.j;
                double zDist = otherAtoms[i].centre.k - atoms[j].centre.k;
                double distance = sqrt(xDist * xDist + yDist * yDist + zDist * zDist);

                energy += (otherAtoms[i].charge * atoms[j].charge) / distance;
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
            int otherNumAtoms = p < leftOver ? atomsPerPlace+1 : atomsPerPlace;
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
    
    if (myRank == 0) {
        printf("energy = %.5f\n", totalEnergy);
    }

    return energy;
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
