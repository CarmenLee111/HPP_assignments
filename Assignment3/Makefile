CFLAGS=-Wall
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: main.o graphics.o
	gcc -O2 -o galsim main.o graphics.o $(LDFLAGS)

main.o: graphics/graphics.h main.c
	gcc $(CFLAGS) $(INCLUDES) -c main.c

graphics.o: graphics/graphics.c graphics/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics/graphics.c

clean:
	rm -f *.o galsim
