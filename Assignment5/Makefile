CFLAGS=-Wall -g 
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: main.c
	gcc $(CFLAGS) -march=native -ffast-math -o galsim main.c -lm


clean:
	rm -rf *.o galsim result.gal testResult.txt galsim.dSYM
	
