. ~/x10-trunk/x10.profile
export TESTFILE=test/3D-tetrapep-cc-pVQZ.inp
for iter in {1..10}
do
  echo "X10_NTHREADS=$X10_NTHREADS bin/pumjarasaayani $TESTFILE"
  bin/pumjarasaayani $TESTFILE
done
