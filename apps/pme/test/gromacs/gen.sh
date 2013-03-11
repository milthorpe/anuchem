cp topol.top box12k.top
genbox -cs spc216.gro -box 5.0 -o box12k.gro -p box12k.top
grompp -v -c box12k.gro -p box12k.top -o box12k.tpr -po box12k.mdp
cp topol.top box27k.top
genbox -cs spc216.gro -box 6.5 -o box27k.gro -p box27k.top
grompp -v -c box27k.gro -p box27k.top -o box27k.tpr -po box27k.mdp
cp topol.top box51k.top
genbox -cs spc216.gro -box 8.0 -o box51k.gro -p box51k.top
grompp -v -c box51k.gro -p box51k.top -o box51k.tpr -po box51k.mdp
cp topol.top box100k.top
genbox -cs spc216.gro -box 10.0 -o box100k.gro -p box100k.top
grompp -v -c box100k.gro -p box100k.top -o box100k.tpr -po box100k.mdp

