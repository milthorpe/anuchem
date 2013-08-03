for iter in {1..10}
do
  echo "X10_NTHREADS=$X10_NTHREADS bin/pumjarasaayani test/water5-cc-pVQZ.inp"
  bin/pumjarasaayani test/water5-cc-pVQZ.inp
done
