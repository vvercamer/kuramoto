#ifndef _RUNGE_H_INCLUDED
#include "toto.h"
#define _RUNGE_H_INCLUDED


double runge4(double x, double *y, double step);	/* Runge-Kutta function */

double f(double x, double y);		/* function for derivative */


#endif
