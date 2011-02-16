/* source : http://www.physics.orst.edu/~rubin/nacphy/ComPhys/DIFFEQ/PRO/srk_c/srk_c.html */

/* Runge Kutta algorithm for first-order differential equations*/

#include <stdio.h>
#include <math.h>
#include "runge.h"
#include "toto.h"

double runge4(double x, double *y, double step, KuramotoStruct* kuramoto)
{
double h=step/2.0, k1, k2, k3, k4;

k1=step*f(x, *y);
k2=step*f(x+h, *y+k1/2.0);
k3=step*f(x+h, *y+k2/2.0);
k4=step*f(x+step, *y+k3);

return(*y+(k1+2*k2+2*k3+k4)/6.0);
}

double  f(double x, double y)
{
return(sin(y));
}


