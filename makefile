CC = g++
CFLAGS = -c -Wall -O3

prog: main.o
	g++ main.o -o prog

main.o: main.cpp
	g++ $(CFLAGS) main.cpp 

clean:
	rm *.o prog
