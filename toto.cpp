#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h> 

#include "runge.h"
#include "gnuplot.h"

int main(int argc, char *argv[])
{
	double tMax = 5;
	double deltaT = 0.01;
	double initial = 1;
	int nbp;
	nbp=(int)floor(tMax/deltaT);
	
	GNUplot gp;
	double* t = (double *) malloc (nbp*sizeof(double));
	double* y = (double *) malloc (nbp*sizeof(double));		

	y[0]=initial;
	algorunge(t,y,nbp, deltaT);
	
	gp.draw(t,y,nbp);

	free(t);
	free(y);	
	
	return 0;

}

