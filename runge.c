/* source : http://www.physics.orst.edu/~rubin/nacphy/ComPhys/DIFFEQ/PRO/srk_c/srk_c.html */

/* Runge Kutta algorithm for first-order differential equations*/

#include <stdio.h>
#include "runge.h"



int algorunge(double * x, double * y, int nbp, double dist)
{
FILE *output = NULL;				/* internal filename */
int j;
 
output=fopen ("runge.dat", "w");	/* external filename */
	if (output != NULL)
	{
		fprintf(output, "0\t%f\n", y[0]);
	
		for (j=1; j<nbp; j++)		/* the time loop */
		{
			x[j]=j*dist;
			y[j]=runge4(x[j], y[j-1], dist);
	
			fprintf(output, "%f\t%f\n", x[j], y[j]);
		}
	
		fclose(output);
	
		return 0;
	}
else 
{
return 1;
}
}

double runge4(double x, double y, double step)
{
double h=step/2.0, k1, k2, k3, k4;

k1=step*f(x, y);
k2=step*f(x+h, y+k1/2.0);
k3=step*f(x+h, y+k2/2.0);
k4=step*f(x+step, y+k3);

return(y+(k1+2*k2+2*k3+k4)/6.0);
}

double  f(double x, double y)
{
return(-y);
}
