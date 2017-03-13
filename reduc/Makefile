CC=mpicxx
CFLAGS = -std=c++11  -O2 #standard c99 , -02 is optimizer
LDFLAGS= 
LNKLIBS= -lm

MPIRUN =mpirun
NPROCS?=1
N?=4
BOOLIC?=1 # 1 for MPI_Allreduce or 0 for the recursive-doubling sum


LIB=libfunction.a
LIBSRC=$(wildcard *.cpp)
LIBOBJ=$(patsubst %.cpp, %.o, $(LIBSRC))
CLEAN_FILES+=$(LIB) $(LIBOBJ)

# Implicit rule for compilation of all objects
%.o:%.cpp
	$(CC) $(CFLAGS) -c $< 


# Rule to build the library
$(LIB): $(LIBOBJ)
	@echo "LIB $@"
	@ar r $@ $^


# Program
all: reduc_mach reduce_zeta


# reduc mach
reduc_mach: reduc_mach.o $(LIB)
	    @$(CC) $(LDFLAGS) -o $@ $^ $(LINKLIBS)
	    $(MPIRUN) -n $(NPROCS) ./$@ $(N) $(BOOLIC)	

# reduc zeta
reduc_zeta: reduc_zeta.o $(LIB)
	    @$(CC) $(LDFLAGS) -o $@ $^ $(LINKLIBS)
	    $(MPIRUN) -n $(NPROCS) ./$@ $(N) $(BOOLIC)	


clean:
	rm -f $(CLEAN_FILES) reduc_zeta reduc_mach
