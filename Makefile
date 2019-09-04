CC=gcc
CXX=g++
CFLAGS=-O0 -g -Wall `pkg-config --cflags cbc` -fsanitize=address
LDFLAGS=-O0 -g -Wall `pkg-config --libs cbc` -fsanitize=address -lm

all:tsp-compact queens queens-lazy tsp-cuts rcpsp

tsp-compact:tsp-compact.o tsp-instance.o
	$(CXX) $(CFLAGS) tsp-compact.o tsp-instance.o -o tsp-compact $(LDFLAGS) -lm

tsp-compact.o:tsp-compact.c tsp-instance.o
	$(CC) $(CFLAGS) -c tsp-compact.c -o tsp-compact.o

tsp-cuts:tsp-cuts.o
	$(CXX) $(CFLAGS) tsp-cuts.o tsp-instance.o spaths.o mincut.o -o tsp-cuts $(LDFLAGS) -lm

tsp-cuts.o:tsp-cuts.c mincut.o tsp-instance.o spaths.o
	$(CC) $(CFLAGS) -c tsp-cuts.c -o tsp-cuts.o

spaths.o:spaths.c spaths.h
	$(CC) $(CFLAGS) -c spaths.c -o spaths.o

mincut.o:mincut.c mincut.h
	$(CC) $(CFLAGS) -c mincut.c -o mincut.o

tsp-instance.o:tsp-instance.c tsp-instance.h
	$(CC) $(CFLAGS) -c tsp-instance.c -o tsp-instance.o

queens:queens.o
	$(CXX) $(CFLAGS) queens.o -o queens $(LDFLAGS) -lm

queens.o:queens.c
	$(CC) $(CFLAGS) -c queens.c -o queens.o

queens-lazy:queens-lazy.o
	$(CXX) $(CFLAGS) queens-lazy.o -o queens-lazy $(LDFLAGS) -lm

queens-lazy.o:queens-lazy.c
	$(CC) $(CFLAGS) -c queens-lazy.c -o queens-lazy.o

rcpsp:rcpsp.o
	$(CXX) $(CFLAGS) rcpsp.o -o rcpsp $(LDFLAGS) -lm

rcpsp.o:rcpsp.c
	$(CC) $(CFLAGS) -c rcpsp.c -o rcpsp.o

clean:
	rm -f *.o tsp-compact queens queens-lazy
