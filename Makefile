CC = g++ 
CFLAGS=-pedantic -Wall -Werror -O2
LDFLAGS=

OBJ = toto.o gnuplot.o gnuplot_i.o
 
all:: toto clean

toto: $(OBJ)
	$(CC) -o toto $(OBJ) $(LDGLAGS)

toto.o: toto.cpp toto.h gnuplot.h gnuplot.h
	$(CC) -c toto.cpp $(CFLAGS)

gnuplot.o: gnuplot.cpp gnuplot.h
	$(CC) -c gnuplot.cpp $(CFLAGS)

gnuplot_i.o: gnuplot_i.c gnuplot_i.h
	$(CC) $(CFLAGS) -c gnuplot_i.c

clean::
	@rm *.o

