CC = g++ 
CFLAGS=-pedantic -Wall -Werror -O2 -I/usr/local/include
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = toto.o  gnuplot_i.o
 
all:: toto clean

toto: $(OBJ)
	$(CC) -L/usr/local/lib -o toto $(OBJ) $(LDFLAGS)

toto.o: toto.cpp toto.h gnuplot_i.h
	$(CC) $(CFLAGS) -c toto.cpp 

gnuplot_i.o: gnuplot_i.c gnuplot_i.h
	$(CC) $(CFLAGS) -c gnuplot_i.c

clean::
	@rm *.o

