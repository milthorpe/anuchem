#ifndef __ELECTROSTATIC_DIRECT_H
#define __ELECTROSTATIC_DIRECT_H
#include <time.h>
#include <iostream>
#include "Atom.h"

using namespace std;

class ElectrostaticDirect {
    static const double SIZE = 80.0;
    long numAtoms;
    Atom *atoms;
    long microTime();
	double randomNoise();

    public:
		ElectrostaticDirect(long);
		~ElectrostaticDirect();
		void run();
		double getEnergy();
};
#endif // ELECTROSTATIC_DIRECT_H


