#!/usr/bin/env bash
# This file is part of ANUChem.
#
#  This file is licensed to You under the Eclipse Public License (EPL);
#  You may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#      http://www.opensource.org/licenses/eclipse-1.0.php
#
# (C) Copyright Josh Milthorpe 2010.
#
# This script runs a calculation for the "bigwater" system of
#   17,132 SPC water molecules (51,396 atoms)
#   using np=256 X10 places (1 place per core)
#   on partition r000n00-c64i4 
#   of the Blue Gene Watson system.
#
# must first allocate partition r000n00-c64i4 using bgpconsole

mpirun -env BG_COREDUMPONEXIT=1 -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition r000n00-c64i4 -mode VN -np 256 bin/pme test/gromacs/conf.gro
