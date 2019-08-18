CC=gcc
CFLAGS=-O0 -g -Wall `pkg-config --cflags cbc`
LDFLAGS=-O0 -g -Wall `pkg-config --libs cbc`

all:tsp-compact queens queens-lazy

tsp-compact:tsp-compact.o tsp-instance.o
	$(CC) $(CFLAGS) tsp-compact.o tsp-instance.o -o tsp-compact $(LDFLAGS) -lm

tsp-compact.o:tsp-compact.c tsp-instance.o
	$(CC) $(CFLAGS) -c tsp-compact.c -o tsp-compact.o

tsp-instance.o:tsp-instance.c tsp-instance.h
	$(CC) $(CFLAGS) -c tsp-instance.c -o tsp-instance.o

queens:queens.o
	$(CC) $(CFLAGS) queens.o -o queens $(LDFLAGS) -lm


queens.o:queens.c
	$(CC) $(CFLAGS) -c queens.c -o queens.o

queens-lazy:queens-lazy.o
	$(CC) $(CFLAGS) queens-lazy.o -o queens-lazy $(LDFLAGS) -lm

queens-lazy.o:queens-lazy.c
	$(CC) $(CFLAGS) -c queens-lazy.c -o queens-lazy.o

clean:
	rm *.o tsp-compact
