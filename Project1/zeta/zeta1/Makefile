CC=mpicxx
CFLAGS =  -O2 -std=c++11#standard c99 , -02 is optimizer
LDFLAGS= 
LNKLIBS= -lm

MPIRUN =mpirun
NPROCS?=1

EXE=exe
EXESRC=main.cpp
EXEOBJ=$(patsubst %.cpp, %.o, $(EXESRC))
CLEAN_FILES+=$(EXE) $(EXEOBJ)

LIB=libfunction.a
LIBSRC=zeta1.cpp
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
$(EXE):	$(EXEOBJ) $(LIB)
	@echo "EXE $@"
	@$(CC) $(LDFLAGS) -o $@ $^ $(LNKLIBS)

run:	$(EXE)
	$(MPIRUN) -n $(NPROCS) ./$<$(N)

clean:
	rm -f $(CLEAN_FILES)
