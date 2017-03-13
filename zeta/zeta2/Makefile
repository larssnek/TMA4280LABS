CC=mpiCC
CFLAGS=-fopenmp -O2 -std=c++11 #Optimization
TARGET=zeta2
LIB=libutils.a
MAINSRC=main.cpp
NUM_THREADS?=1
MPIRUN=mpirun
NPROCS?=1

SRC=$(wildcard *.cpp)

OBJ=$(patsubst %.cpp,%.o,$(SRC))
MAINOBJ=$(patsubst %.cpp,%.o, $(MAINSRC))

LIBOBJ=$(patsubst $(MAINONJ),,$(OBJ))

%.o:%.cpp
	$(CC) $(CFLAGS) -c $<

$(TARGET): $(LIB) $(MAINOBJ)
	$(CC) -o $(TARGET) $(MAINONJ) $(LIB) -lgomp

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)


run:	$(TARGET)
	OMP_NUM_THREADS=$(NUM_THREADS)	$(MPIRUN) -n $(NPROCS) ./$<

clean:
	rm -f $(TARGET) $(LIB) $(OBJ)

list_libobjs: $(LIB)
	@echo "All the objects contained in the stativ library are printed below: \n"
	nm $(LIB)

