#ifndef __ATOM_H
#define __ATOM_H

#include "Point3d.h"

struct Atom {
	Point3d centre;
	char* symbol;
	void* bonds; // dummy field to fill space
	double charge;
};
#endif // ATOM_H


