CFLAGS=-Wall -O3 -pthread
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: main.c
	gcc $(CFLAGS) -march=native -ffast-math  $(INCLUDES) -o galsim main.c graphics/graphics.c $(LDFLAGS)

clean:
	rm -rf *.o galsim result.gal testResult.txt galsim.dSYM
	
