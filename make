all: prog

prog: main.o rungekutta.o
	g++ main.o rungekutta.o -o prog

main.o: main.cpp
	g++ main.cpp -O3

rungekutta.o: rungekutta.cpp
	g++ rungekutta.cpp

clean:
	rm *.o prog
