CC=gcc
CFLAGS=-O0 -g -Wall `pkg-config --cflags cbc`
LDFLAGS=-O0 -g -Wall `pkg-config --libs cbc`

all:tsp-compact

tsp-compact:tsp-compact.o tsp-instance.o
	$(CC) $(CFLAGS) tsp-compact.o tsp-instance.o -o tsp-compact $(LDFLAGS) -lm

tsp-compact.o:tsp-compact.c tsp-instance.o
	$(CC) $(CFLAGS) -c tsp-compact.c -o tsp-compact.o

tsp-instance.o:tsp-instance.c tsp-instance.h
	$(CC) $(CFLAGS) -c tsp-instance.c -o tsp-instance.o

clean:
	rm *.o tsp-compact
