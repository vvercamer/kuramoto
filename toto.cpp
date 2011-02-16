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

	
int main(int argc, char *argv[])
{
	/*définition des conditions expérimentales*/
	double tMax = 5;
	double deltaT = 0.01;
	int nbpts;
	nbpts=(int)floor(tMax/deltaT);	
	double *t = (double *) malloc (nbpts*sizeof(double));
	srand ( time(NULL) );	

	/*définition des oscillatteurs*/
	
	
	int idx1=0,idx2=0;
	int nboscs=4;
	complex_d rComplex(0,0);
	double rayontemp=0, psitemp=0;
	
	struct KuramotoStruct *kuramoto=NULL;
	kuramoto->omega = (double *) malloc (nboscs*sizeof(double));
	kuramoto->K=fmod( rand(), 1);
	kuramoto->rayon = (double *) malloc (nbpts*sizeof(double));
	kuramoto->psi = (double *) malloc (nbpts*sizeof(double));
	
	double *theta = (double *) malloc (nboscs*sizeof(double));
	

	for(idx1 = 0 ; idx1 < nboscs ; idx1++)
	{
		theta[idx1] = fmod( rand(), 2*M_PI);
		kuramoto->omega[idx1] = fmod( rand(), 10000);
	}

	for(idx2 = 0 ; idx2 < nbpts ; idx2++)
	{
		kuramoto->rayon[idx2]=0;
		kuramoto->psi[idx2]=0;
	}
		
	meanField(theta , &rayontemp, &psitemp, nboscs);	
	kuramoto->rayon[0]=rayontemp;
	kuramoto->psi[0]=psitemp;	
	
	printf("r=%f\npsi=%f\n",rayontemp,psitemp);

	/*autres définitions*/
	GNUplot gp;

	/*corps du programme*/

	for (idx2=1; idx2<nbpts; idx2++)         /* the time loop */
        {
         	t[idx2]=idx2*deltaT;
		for(idx1=0 ; idx1<nboscs ; idx1++)
		{       
                	theta[idx1] = runge4(t[idx2], &theta[idx1], deltaT, &kuramoto);
                }
		meanField(theta , &rayontemp, &psitemp, nboscs);
        	kuramoto->rayon[idx2]=rayontemp;
        	kuramoto->psi[idx2]=psitemp;
	}
	/*plot*/
	gp.draw(t,kuramoto->rayon,nbpts);

	/*libération de la mémoire*/
	free(t);
	free(kuramoto->omega);
	free(kuramoto->rayon);
	free(kuramoto->psi);
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

