#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h> 
#include <time.h>
#include <complex>
typedef std::complex<double> complex_d;

#include "toto.h"
// #include "runge.h"
#include "gnuplot.h"

	
int main(int argc, char *argv[])
{
	/*définition des conditions expérimentales*/
	double tMax = 10000;
	double deltaT = 0.1;
	int nbpts;
	nbpts=(int)floor(tMax/deltaT);	
	double *t = (double *) malloc (nbpts*sizeof(double));
	srand ( time(NULL) );	

	/*définition des oscillatteurs*/
	
	
	int idxOsc=0,idxTime=0;
	int nbosc=15;
	complex_d rComplex(0,0);
	double rayontemp=0, psitemp=0;
	

	double *omega = (double *) malloc (nbosc*sizeof(double));
	double K;
	double OMEGA = 700;
	double sigma = 1;
	K = fmod( rand(), 10000)/10000;
	double *rayon = (double *) malloc (nbpts*sizeof(double));
	double *psi = (double *) malloc (nbpts*sizeof(double));
	
	double *theta = (double *) malloc (nbosc*sizeof(double));
	double k1, k2, k3, k4;

	for(idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
	{
		theta[idxOsc] = fmod( rand(), 2*M_PI);
		omega[idxOsc] = OMEGA+sigma*gaussianRand();
	}

	for(idxTime = 0 ; idxTime < nbpts ; idxTime++)
	{
		rayon[idxTime]=0;
		psi[idxTime]=0;
	}
		
	meanField(theta , &rayontemp, &psitemp, nbosc);	
	rayon[0]=rayontemp;
	psi[0]=psitemp;	
	
	printf("r=%f\npsi=%f\n",rayontemp,psitemp);
	printf("K = %f\n",K);
	for(idxOsc=0;idxOsc<nbosc;idxOsc++)
	{
		printf("omega_%d = %f\n",idxOsc+1,omega[idxOsc]);
	}

	/*
	for(idxOsc=0;idxOsc<nbosc;idxOsc++)
	{
		printf("theta_%di = %f\n",idxOsc+1,theta[idxOsc]);
	}
	*/
	
	/*autres définitions*/
	GNUplot gp;

	/*corps du programme*/

	for (idxTime=1; idxTime<nbpts; idxTime++)         /* the time loop */
        {
         	t[idxTime]=idxTime*deltaT;

		for(idxOsc=0 ; idxOsc<nbosc ; idxOsc++)
		{      
			k1=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]);
			k2=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k1/2.0);
			k3=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k2/2.0);
			k4=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k3);
			theta[idxOsc] = theta[idxOsc] + (k1 + 2 * k2 + 2 * k3 + k4)*deltaT/6.0; 
		}

		meanField(theta , &rayontemp, &psitemp, nbosc);
        	rayon[idxTime]=rayontemp;
        	psi[idxTime]=psitemp;
	}
	/*plot*/
	gp.draw(t,rayon,nbpts);

	for(idxOsc=0;idxOsc<nbosc;idxOsc++)
	{
		printf("theta_%df = %f\n",idxOsc+1,theta[idxOsc]);
	}
	
	printf("r=%f\npsi=%f\n",rayontemp,psitemp);
	printf("K = %f\n",K);
	/*libération de la mémoire*/
	free(t);
	free(omega);
	free(rayon);
	free(psi);
	free(theta);
	return 0;

}

int meanField(double *theta, double *rayon, double *psi, int nbosc)
{
	int idxOsc;
	complex_d rComplex(0,0);
	for(idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
	{
		rComplex += exp(complex_d(0,theta[idxOsc]));
	}
	rComplex /= nbosc;
	
	//printf("%f %f\n",real(rComplex),imag(rComplex));
	*rayon = abs(rComplex);
	*psi = arg(rComplex);
	
	return 0;
}

double kuramoto(double omega, double K, double psi, double rayon, double theta)
{
	return(omega+K*rayon*sin(psi-theta));
}

double gaussianRand(void)
{
	double randNum1, randNum2;
	randNum1=fmod(rand(), 100000)/100000;
	randNum2=fmod(rand(), 100000)/100000;
	return sqrt(-2.0*log(randNum1))*cos(2*M_PI);
}


