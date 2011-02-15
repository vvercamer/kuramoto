#ifndef _GNUPLOT_H_INCLUDED

#include <stdio.h>

#define _GNUPLOT_H_INCLUDED

class GNUplot {

public:

	GNUplot();
	~GNUplot();
	void operator ()(const char * command);

	void draw (	double * x,
                	double * y1,
                	double * y2,
                	int n);
	void draw (	double * x,
                	double * y,
                	int n);

protected:

	FILE * 	gnuplotpipe;
	int	serial;
	int	uid;

};

#endif
 
