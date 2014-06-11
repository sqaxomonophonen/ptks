PKGS=sdl2

CC=gcc

CFLAGS=-fopenmp -Ofast -Wall -std=c99 -m64 $(shell pkg-config $(PKGS) --cflags)
#CFLAGS=-O0 -g -Wall -std=c99 -m64 $(shell pkg-config $(PKGS) --cflags)

LINK = -fopenmp $(shell pkg-config $(PKGS) --libs) -lm
#LINK = $(shell pkg-config $(PKGS) --libs) -lm

all: ptks mat.dat

ptks.o: ptks.c
	$(CC) $(CFLAGS) -c ptks.c

ptks: ptks.o
	$(CC) $(LINK) ptks.o -o ptks

mat.dat: floor128.png lited128.png lite128.png wall128.png
	./genmat.py floor128.png: wall128.png: lited128.png:lite128.png > mat.dat

clean:
	rm -rf ptks mat.dat *.o
