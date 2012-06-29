# default shell
gcc TimeStat.c -lm

echo 'w2 DZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water2-cc-pVDZ.inp &> test/water2-cc-pVDZ.out
./time.sh test/water2-cc-pVDZ.out
echo 'w2 TZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water2-cc-pVTZ.inp &> test/water2-cc-pVTZ.out
./time.sh test/water2-cc-pVTZ.out
echo 'w2 QZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water2-cc-pVQZ.inp &> test/water2-cc-pVQZ.out
./time.sh test/water2-cc-pVQZ.out


echo 'w5 DZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water5-cc-pVDZ.inp &> test/water5-cc-pVDZ.out
./time.sh test/water5-cc-pVDZ.out
echo 'w5 TZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water5-cc-pVTZ.inp &> test/water5-cc-pVTZ.out
./time.sh test/water5-cc-pVTZ.out
echo 'w5 QZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water5-cc-pVQZ.inp &> test/water5-cc-pVQZ.out
./time.sh test/water5-cc-pVQZ.out


echo 'w10 DZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water10-cc-pVDZ.inp &> test/water10-cc-pVDZ.out
./time.sh test/water10-cc-pVDZ.out
echo 'w10 TZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water10-cc-pVTZ.inp &> test/water10-cc-pVTZ.out
./time.sh test/water10-cc-pVTZ.out
echo 'w10 QZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water10-cc-pVQZ.inp &> test/water10-cc-pVQZ.out
./time.sh test/water10-cc-pVQZ.out


echo 'w20 DZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water20-cc-pVDZ.inp &> test/water20-cc-pVDZ.out
./time.sh test/water20-cc-pVDZ.out
echo 'w20 TZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water20-cc-pVTZ.inp &> test/water20-cc-pVTZ.out
./time.sh test/water20-cc-pVTZ.out
echo 'w20 QZ'
X10_NTHREADS=1 bin/pumjarasaayani test/water20-cc-pVQZ.inp &> test/water20-cc-pVQZ.out
./time.sh test/water20-cc-pVQZ.out

