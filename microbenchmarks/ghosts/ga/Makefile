GAINC := /scratch/ga/include
GALIB := /scratch/ga/lib
GASOURCE := /scratch/source/ga
F90FLAGS := -I $(GAINC)  -I $(GASOURCE)/armci/src/include
LDFLAGS := -L$(GALIB) -lga -lblas -Wl,-rpath

ghosts.x:	ghosts.F
	mpif90 -o ghosts.x ghosts.F $(F90FLAGS) $(LDFLAGS)

clean:
	rm -rf ghosts.x

