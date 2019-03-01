CFLAGS=-Wall -O3 -march=native -ffast-math -pthread
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm -Xpreprocessor -fopenmp -lomp 

galsim: main.c
	gcc $(CFLAGS) $(INCLUDES) -o galsim main.c graphics/graphics.c $(LDFLAGS)

clean:
	rm -rf *.o galsim result.gal testResult.txt galsim.dSYM
	
