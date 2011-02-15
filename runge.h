#ifndef _RUNGE1ORDRE_H_INCLUDED
#define _RUNGE1ORDRE_H_INCLUDED


int algorunge(double * x, double *y, int nbp, double dist);

double runge4(double x, double y, double step);	/* Runge-Kutta function */

double f(double x, double y);		/* function for derivative */


#endif
