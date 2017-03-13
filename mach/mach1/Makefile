CC=mpiCC
CFLAGS= -O2 -std=c++11 #Optimization
TARGET=mach1
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
	$(CC) -o $(TARGET) $(MAINONJ) $(LIB)

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)


run:	$(TARGET)
	$(MPIRUN) -n $(NPROCS) ./$<	
#OMP_NUM_THREADS=$(NUM_THREADS)	$(MPIRUN) -n $(NPROCS) ./$<

clean:
	rm -f $(TARGET) $(LIB) $(OBJ)

list_libobjs: $(LIB)
	@echo "All the objects contained in the stativ library are printed below: \n"
	nm $(LIB)



# mach1: main.o mach1.o 
#	mpiCC -o mach1 main.o mach1.o
#
#main.o: main.cpp mach1.h
#	mpiCC -c main.cpp
#
#mach1.o: mach1.cpp mach1.h
#	mpiCC -c mach1.cpp

