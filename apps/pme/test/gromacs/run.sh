#!/bin/bash
part[0]=12;  s[0]=5.0
part[1]=27;  s[1]=6.5
part[2]=51;  s[2]=8.0
part[3]=100; s[3]=10.0
part[4]=194; s[4]=12.5
part[5]=497; s[5]=17.0
part[6]=1067; s[6]=22.0

for i in 0 1 2 3 4 5 6
do
    prefix=box${part[i]}k
    grompp -v -f grompp.mdp -c ${prefix}.gro -p ${prefix}.top -o ${prefix}.tpr -po ${prefix}.mdp
    mdrun -s ${prefix}.tpr -c ${prefix}.gro -nsteps 1000 -g ${prefix}.log -nt 1
done

