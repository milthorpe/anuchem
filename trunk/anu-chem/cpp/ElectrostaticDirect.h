#ifndef __ELECTROSTATIC_DIRECT_H
#define __ELECTROSTATIC_DIRECT_H
#include <time.h>
#include <iostream>
#include "Atom.h"

using namespace std;

static const double SIZE = 80.0;

int numAtoms;
Atom* atoms;

int myRank;
int nTasks;

int atomsPerPlace;
int leftOver;
int myNumAtoms;
int myStart;
int myEnd;

static long microTime();
static double randomNoise();

static double getEnergy(Atom* atoms);

#endif // ELECTROSTATIC_DIRECT_H


