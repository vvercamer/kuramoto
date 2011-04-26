#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h> 
#include <time.h>
#include <complex>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
typedef std::complex<double> complex_d;

#include "toto.h"
#include "gnuplot_i.h"

	
int main(int argc, char *argv[])
{
	double deltaT = 0.1;
	int nbsamples=10000;
	int nbosc = 100;
	int nbK=10;	
	
	if (argc == 5)
	{
		deltaT = strtod(argv[1],NULL);
		nbsamples = (int)strtol(argv[2], NULL, 10);
		nbosc = (int)strtol(argv[3], NULL, 10);
		nbK = (int)strtol(argv[4], NULL, 10);
	}
	else if (argc == 1)
	{
	}
	else if (argc == 2)
	{
		if (strcmp(argv[1],"-h\n"))
		{
			printf("utilisation :\n ./toto deltaT nbsamples nbosc nbK\n");
		}
		return 0;
	}
	else 
	{
		printf("problème sur les paramètres d'entrée\n");
		return 1;
	}
	/*lecture des paramètres */
/*	if (argc>1)
	{
		int idx;
		for(idx = 1 ; idx<argc ; idx++)
		{
			 *argv[idx];
			switch(s)
			{
			case "deltaT":
			  deltaT = strtod(agrv[idx+1],NULL);
			  break;
			case "nbsamples":
			  nbsamples = (int)strtol(argv[idx+1], NULL, 10);
			  break;
			case "nbosc":
			  nbosc = (int)strtol(argv[idx+1], NULL, 10);
			  break;
			case "nbK":
			  nbK = (int)strtol(argv[idx+1], NULL, 10);
			  break;
			default:
			  break;
			}
		}
	}
*/
	/*définition des conditions expérimentales*/
	double *temps = (double *) malloc (nbsamples*sizeof(double));
	srand ( time(NULL) );	

	/*définition des oscillatteurs*/
	
	int idxOsc, idxTime, idxK;
	
	double rayontemp=0, psitemp=0;
	
	double *omega = (double *) malloc (nbosc*sizeof(double));
	double *theta = (double *) malloc (nbosc*sizeof(double));
	double K = 0;
	double OMEGA = 0;
	//double sigma = 0.3;
	double *rayon = (double *) malloc (nbsamples*sizeof(double));
	double *psi = (double *) malloc (nbsamples*sizeof(double));
	double *rayonmoyen = (double *) malloc (nbsamples*sizeof(double));

	double *rayonstable = (double *) malloc (nbK*sizeof(double));
	double *invrayonstable = (double *) malloc (nbK*sizeof(double));
	double *Kvect = (double *) malloc (nbK*sizeof(double));
	
	double k1, k2, k3, k4;
	
	const gsl_rng_type * randType;
		gsl_rng * r;
		gsl_rng_env_setup();
		randType = gsl_rng_default;
		r = gsl_rng_alloc (randType);
	
	for(idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
	{
		//omega[idxOsc] = OMEGA + gsl_ran_gaussian(r,sigma);
		omega[idxOsc] = OMEGA + gsl_ran_cauchy(r,0.25);
	}
	

	for(idxK = 0 ; idxK < nbK ; idxK++)
	{
		rayonstable[idxK]=0;
        int idxTimeStart;
        idxTimeStart=(int)floor(3*nbsamples/4.);

		if (nbK == 1)
		{
			printf("entrer la valeur de K\n");
			scanf("%lf",&K);
		}
		else
		{
			K = (double)idxK*2/(nbK-1);
		}
	
		for(idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
		{
			theta[idxOsc] = fmod( rand(), 2*M_PI)-M_PI;
		}
	
		for(idxTime = 0 ; idxTime < nbsamples ; idxTime++)
		{
			rayon[idxTime]=0;
			psi[idxTime]=0;
		}
			
		meanField(theta , &rayontemp, &psitemp, nbosc);	
		rayon[0]=rayontemp;
		psi[0]=psitemp;	
		
		printf("K = %f\n",K);
		
		
		/*corps du programme*/
	
		for (idxTime = 1 ; idxTime < nbsamples ; idxTime++)         /* the time loop */
	    {
	      	temps[idxTime]=idxTime*deltaT;
	
			for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
			{      
				/*méthode de runge-kutta d'ordre 4*/
				k1=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]);
				k2=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k1/2.0);
				k3=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k2/2.0);
				k4=deltaT*kuramoto(omega[idxOsc], K, psi[idxTime-1], rayon[idxTime-1],theta[idxOsc]+deltaT*k3);
				theta[idxOsc] = theta[idxOsc] + (k1 + 2 * k2 + 2 * k3 + k4)*deltaT/6.0; 
			}
			
			meanField(theta , &rayontemp, &psitemp, nbosc);
        	rayon[idxTime]=rayontemp;
        	psi[idxTime]=psitemp;
			if (idxTime > idxTimeStart)
			{
				rayonstable[idxK]+=rayontemp;
			}		

		}
		
		rayonstable[idxK]/=nbsamples-idxTimeStart+1;
		invrayonstable[idxK]= 1/(1 - rayonstable[idxK]*rayonstable[idxK]);//pour formule 4.7 
		Kvect[idxK]=K;	
	}	


	if (nbK == 1)
	{
		rayonmoyen[0]=rayon[0];
		for (idxTime = 1 ; idxTime < nbsamples ; idxTime++)
		{
			rayonmoyen[idxTime] = rayonmoyen[idxTime-1] + rayon[idxTime];
		}
		for (idxTime = 0 ; idxTime < nbsamples ; idxTime++)
        {
            rayonmoyen[idxTime]/=idxTime+1;
        }
	}


	/*définitions pour le graphique*/
	gnuplot_ctrl * gp1;
	gnuplot_ctrl * gp2;
	gp1 = gnuplot_init();
	gp2 = gnuplot_init();

	char titre[256];
#if defined ( __APPLE__ )
    gnuplot_cmd(gp1, "set terminal x11 persist");
	gnuplot_cmd(gp2, "set terminal x11 persist");
#else
	gnuplot_cmd(gp1, "set terminal wxt persist");
	gnuplot_cmd(gp2, "set terminal wxt persist");
#endif
	gnuplot_setstyle(gp1, "lines");	
	gnuplot_set_ylabel(gp1, "r");
	gnuplot_setstyle(gp2, "lines");
    gnuplot_set_ylabel(gp2, "r");

	if (nbK == 1)
	{
		gnuplot_set_xlabel(gp1, "t");
		gnuplot_cmd(gp1, "set yrange [-0.05:1.05]");
		sprintf(titre,"evolution de r(t) pour K = %f", K);
		gnuplot_plot_xy(gp1, temps, rayon, nbsamples, titre) ;
		sprintf(titre,"evolution de rmoyen(t) pour K = %f", K);
		gnuplot_plot_xy(gp1, temps, rayonmoyen, nbsamples, titre) ;
	}
	else
	{
		gnuplot_set_xlabel(gp1, "K");
        gnuplot_cmd(gp1, "set yrange [-0.05:1.05]");
		gnuplot_plot_xy(gp1, Kvect, rayonstable, nbK, "evolution de rstable(K)") ;
	
		gnuplot_set_xlabel(gp2, "K");
        gnuplot_cmd(gp2, "set yrange [-0.05:10.05]");
        gnuplot_plot_xy(gp2, Kvect, invrayonstable, nbK, "evolution de 1/(1-stable(K)^2)") ;
	}

	printf("r=%f\npsi=%f\n",rayontemp,psitemp);
	
	/*libération de la mémoire*/
	free(temps);
	free(omega);
	free(rayon);
	free(psi);
	free(theta);
	gnuplot_close(gp1);
	gnuplot_close(gp2);
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
	
	*rayon = abs(rComplex);
	*psi = arg(rComplex);
	
	return 0;
}

double kuramoto(double omega, double K, double psi, double rayon, double theta)
{
	return(omega+K*rayon*sin(psi-theta));
}


