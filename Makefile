


CC = g++
CFLAGS  = -O3

all: mainPar mainSer

mainPar:  mainPar.cpp 
	mpic++ mainPar.cpp $(CFLAGS) -o mainPar.o 


mainSer:  mainSer.cpp
	$(CC) mainSer.cpp $(CFLAGS) -o mainSer.o 

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~
