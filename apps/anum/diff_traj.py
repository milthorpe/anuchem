# IPython log file

import struct
import numpy as np
import math
import sys

vector3d = np.dtype([('x', '>d'), ('y', '>d'), ('z', '>d')])
dt = np.dtype([('index', '>i4'), ('species', '>i4'), ('pos', vector3d), ('v', vector3d)])

if (len(sys.argv) < 5):
    sys.exit('usage: diffs numSnapshots timestepsPerSnapshot dir1 dir2')
numSnapshots=int(sys.argv[1])
ts=int(sys.argv[2])
dir1=sys.argv[3]
dir2=sys.argv[4]

print('differences for ' + str(numSnapshots) + ' snapshots between ' + dir1 + ' and ' + dir2)

difffile = open('diffs.txt','w+')
difffile.write('#time  rms(position) rms(velocity)\n');

for snap in range(1,numSnapshots+1):
    time=snap*ts
    print(str(time))

    snapfile='snapshot_'+str(time)+'.dat'
    file1=dir1+'/'+snapfile
    file2=dir2+'/'+snapfile

    fo1 = open(file1, 'r+')
    fo2 = open(file2, 'r+')

    time1 = struct.unpack('>i', fo1.read(4))
    if (time1[0] != time):
        sys.exit('found invalid snapshot ' + str(file1) + 'for time ' + str(time1) + ' expected ' + str(time))
    time2 = struct.unpack('>i', fo2.read(4))
    if (time2[0] != time):
        sys.exit('found invalid snapshot ' + str(file2) + 'for time ' + str(time2) + ' expected ' + str(time))
    natoms1 = struct.unpack('>i', fo1.read(4))
    natoms2 = struct.unpack('>i', fo2.read(4))


    atoms1 = np.fromfile(fo1, dtype=dt, count=natoms1[0])
    atoms1 = sorted(atoms1, key=lambda dt:dt[0])

    atoms2 = np.fromfile(fo2, dtype=dt, count=natoms2[0])
    atoms2 = sorted(atoms2, key=lambda dt:dt[0])

    msex = 0.0
    msev = 0.0
    normx = 0.0
    normv = 0.0
    j=0
    for i in range(0,natoms1[0]):
        while (atoms2[j][0] < atoms1[i][0]):
            j = j+1
        if (j >= natoms2[0]):
            break
        while (atoms1[i][0] < atoms2[j][0]):
            i = i+1  
        magx = atoms2[j][2][0]*atoms2[j][2][0]+atoms2[j][2][1]*atoms2[j][2][1]+atoms2[j][2][2]*atoms2[j][2][2]
        normx += magx
        diffx = atoms2[j][2][0]-atoms1[i][2][0]
        diffy = atoms2[j][2][1]-atoms1[i][2][1]
        diffz = atoms2[j][2][2]-atoms1[i][2][2]
        diff = diffx*diffx+diffy*diffy+diffz*diffz
        msex += diff
        magv = atoms2[j][3][0]*atoms2[j][3][0]+atoms2[j][3][1]*atoms2[j][3][1]+atoms2[j][3][2]*atoms2[j][3][2]
        normv += magv
        diffvx = atoms2[j][3][0]-atoms1[i][3][0]
        diffvy = atoms2[j][3][1]-atoms1[i][3][1]
        diffvz = atoms2[j][3][2]-atoms1[i][3][2]
        diffv = diffvx*diffvx+diffvy*diffvy+diffvz*diffvz
        msev += diffv

    rmsx = math.sqrt(msex/normx)
    rmsv = math.sqrt(msev/normv)
    difffile.write('{} {} {}\n'.format(time,rmsx,rmsv))

    fo1.close()
    fo2.close()

difffile.close()


