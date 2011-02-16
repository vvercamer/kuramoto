#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h> 
#include <time.h>
#include <complex>
typedef std::complex<double> complex_d;

#include "toto.h"
#include "runge.h"
#include "gnuplot.h"

int meanField(double *theta, double *rayon, double *psi, int nboscs);

int main(int argc, char *argv[])
{
	/*définition des conditions expérimentales*/
	double tMax = 5;
	double deltaT = 0.01;
	int nbpts;
	nbpts=(int)floor(tMax/deltaT);	
	double *t = (double *) malloc (nbpts*sizeof(double));
	double *y = (double *) malloc (nbpts*sizeof(double));		
	srand ( time(NULL) );	

	/*définition des oscillatteurs*/
struct Kuramoto{

};
	int idx1=0,idx2=0;
	int nboscs=4;
	complex_d rComplex(0,0);
	double rayon=0, psi=0;
	double *theta = (double *) malloc (nboscs*sizeof(double));

	for(idx1=0 ; idx1<nboscs ; idx1++)
	{
		theta[idx1] = fmod( rand() , 2*M_PI);
	}
		
	meanField(theta , &rayon, &psi, nboscs);	
	
	printf("r=%f\npsi=%f\n",rayon,psi);

	/*autres définitions*/
	GNUplot gp;

	/*corps du programme*/

	for(idx1=0 ; idx1<nboscs ; idx1++)
	{
		y[0] = theta[idx1];
		printf("theta i initial = %f\n",y[0]);
		//algorunge(t,y,nbpts, deltaT);
		for (idx2=1; idx2<nbpts; idx2++)         /* the time loop */
                {
                        t[idx2]=idx2*deltaT;
                        y[idx2]=runge4(t[idx2], y[idx2-1], deltaT);
                }

	}
	/*plot*/
//	gp.draw(t,y,nbpts);

	/*libération de la mémoire*/
	free(t);
	free(y);	
	free(theta);
	return 0;

}

int meanField(double *theta, double *rayon, double *psi, int nboscs)
{
	int idx1;
	complex_d rComplex(0,0);
	for(idx1 = 0 ; idx1 < nboscs ; idx1++)
	{
		rComplex += exp(complex_d(0,theta[idx1]));
	}
	rComplex /= nboscs;
	*rayon = abs(rComplex);
	*psi = arg(rComplex);
	
	return 0;
}

