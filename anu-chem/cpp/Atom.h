#ifndef __ATOM_H
#define __ATOM_H

#include "Point3d.h"

struct Atom {
	Point3d centre;
	double charge;
    Point3d force;
};
#endif // ATOM_H


