#include "ElectrostaticDirect.h"
#include <stdlib.h>
#include <sys/time.h> 
#include <math.h>
#include <stdio.h>

ElectrostaticDirect::ElectrostaticDirect(long N) {
	numAtoms = N;
	atoms = new Atom[numAtoms];
}

ElectrostaticDirect::~ElectrostaticDirect() {
	delete atoms;
}

long ElectrostaticDirect::microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}

void ElectrostaticDirect::run() {
    int gridSize = (int)(ceil(cbrt(numAtoms)));
	cout << gridSize << '\n';
    int gridPoint = 0; // running total of assigned grid points
	for (int i=0; i<numAtoms; i++) {
		int gridX = gridPoint / (gridSize * gridSize);
        int gridY = (gridPoint - (gridX * gridSize * gridSize)) / gridSize;
        int gridZ = gridPoint - (gridX * gridSize * gridSize) - (gridY * gridSize);
        double x = (gridX + 0.5 + randomNoise()) * (SIZE / gridSize);
        double y = (gridY + 0.5 + randomNoise()) * (SIZE / gridSize);
        double z = (gridZ + 0.5 + randomNoise()) * (SIZE / gridSize);
        atoms[i].centre.i = x;
		atoms[i].centre.j = y;
		atoms[i].centre.k = z;
		atoms[i].charge = i%2==0?1.0:-1.0;
        gridPoint++;
	}

    long start = microTime();
	getEnergy();
    long stop = microTime();
    printf("%.3g s\n", (stop-start) / 1.0e6);
}

double ElectrostaticDirect::randomNoise() {
	return 0.0;
}

double ElectrostaticDirect::getEnergy() {
    double energy = 0.0;

    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < i; j++) {

            double xDist = atoms[i].centre.i -atoms[j].centre.i;
            double yDist = atoms[i].centre.j -atoms[j].centre.j;
            double zDist = atoms[i].centre.k -atoms[j].centre.k;

            double distance = sqrt(xDist * xDist +
                    yDist * yDist +
                    zDist * zDist);

            energy += (atoms[i].charge * atoms[j].charge) / distance;
        }
    }
    printf("energy = %.5f\n", energy);

    return energy;
}

int main(int argc, char* argv[])
{
	long N = 30000;
	if (argc > 1) {
		char* end;
		N = strtol(argv[1], &end, 10);
	}
    ElectrostaticDirect* electro = new ElectrostaticDirect(N);
    electro->run();
    return 0;
}
