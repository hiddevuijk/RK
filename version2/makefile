CC = g++
CFLAGS = -c -Wall -O3

outfile = test

all: $(outfile)

$(outfile): main.o integrator.o eulerInt.o rkInt.o
	$(CC)  main.o integrator.o eulerInt.o rkInt.o -o $(outfile)

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp 

integrator.o: integrator.cpp
	$(CC) $(CFLAGS) integrator.cpp

eulerInt.o: eulerInt.cpp
	$(CC) $(CFLAGS) eulerInt.cpp

rkInt.o: rkInt.cpp
	$(CC) $(CFLAGS) rkInt.cpp

clean:
	rm *.o $(outfile)
