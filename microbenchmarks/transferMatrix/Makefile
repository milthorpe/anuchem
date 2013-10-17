X10CPP := x10c++
X10CPPFLAGS := -x10rt mpi -O -NO_CHECKS #-gpt

run:    transferMatrix
	CPUPROFILE=my.prof mpiexec -n 2 -x X10_NTHREADS=1 ./transferMatrix

transferMatrix:	BenchmarkTransferMatrix.x10
	$(X10CPP) $(X10CPPFLAGS) -x10lib $(X10_HOME)/../x10.gml/native_gml.properties BenchmarkTransferMatrix.x10 -o transferMatrix -post '# # # -L/home/josh/GotoBLAS2 -lgoto2 -lgfortran'

clean:
	rm -rf transferMatrix BenchmarkTransferMatrix.cc BenchmarkTransferMatrix.h *.dat *.ps harness

