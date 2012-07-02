sed -n '/rawtime/,/rawtime/p' $1  > temp
sed -e 's/^.*rawtime//' -e 's!rawtime.*!!' temp > $1.time
rm temp
./a.out $1.time

