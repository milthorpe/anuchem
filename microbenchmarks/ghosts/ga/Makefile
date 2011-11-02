F90FLAGS := -I /scratch/ga/include  -I /scratch/source/ga/armci/src/include
LDFLAGS := -L/scratch/ga/lib -lga -lblas -Wl,-rpath

ghosts.x:	ghosts.F
	mpif90 -o ghosts.x ghosts.F $(F90FLAGS) $(LDFLAGS)

clean:
	rm -rf ghosts.x

