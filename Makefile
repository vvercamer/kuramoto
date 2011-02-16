CC = g++ 
CFLAGS=-pedantic -Wall -Werror -O
LDFLAGS=

OBJ = toto.o runge.o gnuplot.o 
 
all:: toto clean

toto: $(OBJ)
	$(CC) -o toto $(OBJ) $(LDGLAGS)

toto.o: toto.cpp toto.h runge.h gnuplot.h
	$(CC) -c toto.cpp $(CFLAGS)

runge.o: runge.c
	$(CC) -c runge.c $(CFLAGS)
	
gnuplot.o: gnuplot.cpp gnuplot.h
	$(CC) -c gnuplot.cpp $(CFLAGS)

clean::
	@rm *.o

