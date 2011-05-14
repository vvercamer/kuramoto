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
#include <gsl/gsl_fit.h>
#include <gsl/gsl_histogram.h>
typedef std::complex<double> complex_d;

#include "toto.h"
#include "gnuplot_i.h"


int main(int argc, char *argv[])
{
	double deltaT = 0.1;
	int nbsamples = 10000;
	int nbosc = 100;
	int nbrand = 5;
	int nbK = 10;

	if (argc == 6) {
		deltaT = strtod(argv[1],NULL);
		nbsamples = (int)strtol(argv[2], NULL, 10);
		nbosc = (int)strtol(argv[3], NULL, 10);
		nbrand = (int)strtol(argv[4], NULL, 10);
		nbK = (int)strtol(argv[5], NULL, 10);
	}
	else if (argc == 1) {
	}
	else if (argc == 2) {
		if (strcmp(argv[1],"-h\n")) {
		printf("utilisation :\n ./toto deltaT nbsamples nbosc nbrand nbK\n");
		}
		return 0;
	}
	else {
		printf("problème sur les paramètres d'entrée\n");
		return 1;
	}


	/*Déclaration des variables*/

	/*Définition des conditions expérimentales*/
	double *temps = (double *) malloc (nbsamples*sizeof(double));
	srand(time(NULL));

	/*Définition des oscillatteurs*/
	int idxOsc, idxTime, idxK, idxRand;

	double rayontemp = 0, psitemp = 0;

	double *omega = (double *) malloc (nbosc*sizeof(double));
	double *theta = (double *) malloc (nbosc*sizeof(double));
	double K = 0;
	double Kmax = 2;
	double OMEGA = 0;
	double sigma = 0.2;
	double subcrit = 0;
	double *rayon = (double *) malloc (nbsamples*sizeof(double));
	double *psi = (double *) malloc (nbsamples*sizeof(double));

	double *rayonmoyen = (double *) malloc (nbsamples*sizeof(double));
	double *rayonmoyenRand = (double *) malloc (nbsamples*sizeof(double)); /*moyenne sur les réalisations*/
	double *rayonstable = (double *) malloc (nbK*sizeof(double));
//	double *invrayonstable = (double *) malloc (nbK*sizeof(double));
	double *Kvect = (double *) malloc (nbK*sizeof(double));

	/*Détermination de la distribution des pulsations propres*/
	int nbw = 1001;
	int idxw;
	double *w = (double *) malloc (nbw*sizeof(double));
	double *N = (double *) malloc (nbw*sizeof(double));

	/*Détermination de Kc et de beta*/
	double Kc, beta;
	Kc = 2 / (M_PI * gsl_ran_cauchy_pdf(subcrit, sigma));

	printf("Kc = %f\n", Kc);
	beta = pow(Kc, 3) * (1 - 4 * pow(subcrit, 2) / (pow(sigma, 2) * (1 + pow(subcrit / sigma, 2))))/ (8 * pow(sigma, 3) * pow(1 + pow(subcrit / sigma, 2), 2));
	printf("beta = %f\n", beta);

	/*Déclaration des variables nécessaires à la détermination du rayon asymptotique rayonInfini*/
	int idxC;
	double rayonMax = 0;
	double *Tc = (double *) malloc (nbK*sizeof(double));
	double *rayonCut = (double *) malloc (nbsamples*sizeof(double));
	double *rayoninfini = (double *) malloc (nbK*sizeof(double));


	/*Déclaration des variables nécessaires pour Runge-Kutta 4*/

	double k1, k2, k3, k4;


	/*Définition des paramètres pour les histogrammes*/
	double wmax = 0.8;
	double thetamax = M_PI + 0.2;
	size_t nhist = 100;
	int nbhist = (int) nhist;

	double *hist = (double *) malloc (nbhist*sizeof(double));
	double *histthetadeb = (double *) malloc (nbhist*sizeof(double));
	double *histthetafin = (double *) malloc (nbhist*sizeof(double));

	gsl_histogram *htheta = gsl_histogram_alloc(nhist);
	gsl_histogram_set_ranges_uniform (htheta, -thetamax, thetamax);
		


	/*Boucle sur les valeurs de K*/
	for (idxK = 0 ; idxK < nbK ; idxK++) {

		/*Détermination de K et de Kvect*/
		if (nbK == 1) {
			printf("entrer la valeur de K\n");
			scanf("%lf", &K);
		}
		else {
			K = (double)idxK * Kmax / (nbK - 1);
		}

		Kvect[idxK] = K;
		printf("K = %f\n",K);

		/*Détermination des paramètres de calcul de rstable*/
		rayonstable[idxK] = 0;
		int idxTimeStart;
		idxTimeStart = (int)floor(3 * nbsamples / 4.);

		/*Initialisation de rayonmoyenRand*/
		for (idxTime=1 ; idxTime < nbsamples ; idxTime++) {
			rayonmoyenRand[idxTime] = 0;
		}


		/*Boucle sur les réalisations*/
		for (idxRand = 0 ; idxRand < nbrand ; idxRand++) {
		
			/*Détermination des pulsations propres des oscillateurs*/
			gsl_rng * r;
			const gsl_rng_type * T;
			gsl_rng_env_setup();
			gsl_rng_default_seed = rand();
			T = gsl_rng_default;
			r = gsl_rng_alloc (T);
//			printf ("seed = %lu\n", gsl_rng_default_seed);

			for (idxOsc = 0 ; idxOsc < (nbosc / 2) ; idxOsc++) {
//				omega[idxOsc] = OMEGA + gsl_ran_gaussian(r,sigma) + subcrit;
				omega[idxOsc] = OMEGA + gsl_ran_cauchy(r,sigma) + subcrit;
			}
		
			for (idxOsc = (nbosc / 2) ; idxOsc < nbosc ; idxOsc++) {
//				omega[idxOsc] = OMEGA + gsl_ran_gaussian(r,sigma) - subcrit;
				omega[idxOsc] = OMEGA + gsl_ran_cauchy(r,sigma) - subcrit;
			}

			/*Histogramme des pulsations*/
			if (idxRand == 0) {
				gsl_histogram *h = gsl_histogram_alloc(nhist);
				gsl_histogram_set_ranges_uniform (h, -wmax, wmax);
		
				for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
					gsl_histogram_increment (h, omega[idxOsc]);
				}

				for (idxOsc = 0 ; idxOsc < nbhist ; idxOsc++) {
					hist[idxOsc] = gsl_histogram_get (h, idxOsc);
				}
			}

			/*Initialisation de theta, de rayon et de psi*/
			for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
				theta[idxOsc] = fmod(rand(), 2 * M_PI) - M_PI;
			}

			for (idxTime = 0 ; idxTime < nbsamples ; idxTime++) {
				rayon[idxTime] = 0;
				psi[idxTime] = 0;
			}

			meanField(theta , &rayontemp, &psitemp, nbosc);
			rayon[0] = rayontemp;
			psi[0] = psitemp;


			/*Corps du programme*/

			for (idxTime = 1 ; idxTime < nbsamples ; idxTime++) {
				temps[idxTime] = idxTime*deltaT;

				for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
					
					/*Détermination des theta par la méthode de runge-kutta d'ordre 4*/
					k1 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc]);
					k2 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k1 / 2.0);
					k3 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k2 / 2.0);
					k4 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k3);
					theta[idxOsc] = fmod(theta[idxOsc] + (k1 + 2 * k2 + 2 * k3 + k4) * deltaT / 6.0 + M_PI, 2*M_PI) - M_PI;
				}

				/*Histogramme des valeurs initiales des theta*/
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == 1) {
					for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
						gsl_histogram_increment (htheta, theta[idxOsc]);
					}

					for (idxOsc = 0 ; idxOsc < nbhist ; idxOsc++) {
						histthetadeb[idxOsc] = gsl_histogram_get (htheta, idxOsc);
					}
				}

				/*Histogramme des valeurs finales des theta*/
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == nbsamples - 1) {
					for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
						gsl_histogram_increment (htheta, theta[idxOsc]);
					}

					for (idxOsc = 0 ; idxOsc < nbhist ; idxOsc++) {
						histthetafin[idxOsc] = gsl_histogram_get (htheta, idxOsc);
					}
				}
	
				/*Détermination de r et psi par l'approche champ moyen*/
				meanField(theta , &rayontemp, &psitemp, nbosc);
				rayon[idxTime] = rayontemp;
				psi[idxTime] = fmod(psitemp, 2 * M_PI);

				/*Détermination de rstable et rayonmoyenRand*/
				if (idxTime > idxTimeStart)
					rayonstable[idxK] += rayontemp;

				rayonmoyenRand[idxTime] += rayontemp;
			}
		}

		/*Fin de la boucle sur les réalisations*/

		/*Détermination de rayonmoyenRand, moyenne de r sur plusieurs réalisations*/
		for (idxTime = 1 ; idxTime < nbsamples ; idxTime++) {
			rayonmoyenRand[idxTime]/=nbrand;
		}

		/*Détermination de rayonstable et de son inverse*/
		rayonstable[idxK] /= (nbsamples - idxTimeStart + 1) * nbrand;
//		invrayonstable[idxK] = 1 / (1 - rayonstable[idxK] * rayonstable[idxK]);	//pour formule 4.7

		/*Définition de rayonMax, le maximum de r*/
		for (idxTime = 1 ; idxTime < nbsamples ; idxTime++) {
			if (rayonmoyenRand[idxTime] > rayonMax)
				rayonMax = rayonmoyenRand[idxTime];
		}

		/*Définition du temps critique comme étant le temps pour lequel on atteint 90% de la valeur maximale*/
		idxTime = 0;
		while (rayonmoyenRand[idxTime] < 0.9 * rayonMax) 
			idxTime++;

		idxC = idxTime;
		Tc[idxK] = idxC*deltaT;

		/*Détermination de rayonInfini*/
		if (idxC < nbsamples) {
			for (idxTime = idxC; idxTime < nbsamples; idxTime++) {
				rayonCut[idxTime - idxC] = rayonmoyenRand[idxTime];
			}
			rayoninfini[idxK] = gsl_stats_mean(rayonCut, 1, nbsamples - idxC + 1);
//			printf("%g pour K = %d \n",rayoninfini[idxK], idxK);
		}
		else {
//			printf("%g pour K= %d \n", rayoninfini[idxK], idxK);
			rayoninfini[idxK] = gsl_stats_mean(rayonmoyenRand, 1, nbsamples);
			/*car sinon on obtiendrait des valeurs nulles de rayonInfini pour certaines valeurs de K < Kc */
		}

		/*Détermination de rayonmoyen, comme moyenne temporelle du signal entier*/
		if (nbK == 1) {
			rayonmoyen[0] = rayon[0];
			for (idxTime = 1 ; idxTime < nbsamples ; idxTime++) {
				rayonmoyen[idxTime] = rayonmoyen[idxTime - 1] + rayon[idxTime];
			}
			for (idxTime = 0 ; idxTime < nbsamples ; idxTime++) {
				rayonmoyen[idxTime] /= idxTime + 1;
			}
		}
	}


	/*Définition pour les logarithmes*/
	if (nbK > 4) {
		logarithme(Kc, Tc, nbK, Kmax, Kvect, rayonstable);
	}

	/*Détermination de la distribution des pulsations propres*/
	for (idxw=0 ; idxw < nbw ; idxw++) {
		w[idxw] = (-(nbw - 1) / 2 + idxw) * wmax * 2 / nbw;
		N[idxw] = gsl_ran_cauchy_pdf(w[idxw] - subcrit, sigma) / 2 + gsl_ran_cauchy_pdf(w[idxw] + subcrit, sigma) / 2;
	}

	/*Détermination des vecteurs des abscisses pour les distributions et les histogrammes*/
	double *Nosc = (double *) malloc (nbosc*sizeof(double));
	double *Nhist = (double *) malloc (nbhist*sizeof(double));
	double *Nhisttheta = (double *) malloc (nbhist*sizeof(double));

	for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
		Nosc[idxOsc] = idxOsc;
	}
	for (idxOsc = 0 ; idxOsc < nbhist ; idxOsc++) {
		Nhist[idxOsc] = idxOsc * 2 * wmax / nbhist - wmax;
		Nhisttheta[idxOsc] = idxOsc * 2 * thetamax / nbhist - thetamax;	
	}



	/*Définitions pour le graphique*/
	gnuplot_ctrl * gp;
	char titre[256];



	/*Tracé de la distribution théorique des pulsations*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 0 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 0 persist");
#endif
//	gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//	gnuplot_cmd(gp, "set output 'distribution.jpeg'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "pulsation w");
	gnuplot_set_ylabel(gp, "distribution g(w)");
	gnuplot_plot_xy(gp, w, N, nbw,"distribution des pulsations propres des oscillateurs");



	/*Tracé de l'histogramme des valeurs des pulsations*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 1 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 1 persist");
#endif
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set yrange [0:%d]", (int)nbosc / 10);
	gnuplot_set_xlabel(gp, "pulsation w");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	gnuplot_plot_xy(gp, Nhist, hist, nbhist,"histogramme des valeurs des pulsations");	



	/*Tracé de l'histogramme des valeurs initiales des theta*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 2 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 2 persist");
#endif
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set yrange [0:%d]", (int)nbosc / 2);
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	gnuplot_set_xlabel(gp, "phase theta");
	gnuplot_plot_xy(gp, Nhisttheta, histthetadeb, nbhist,"histogramme des valeurs initiales des theta");	



	/*Tracé de l'histogramme des valeurs finales des theta*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 3 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 3 persist");
#endif
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set yrange [0:%d]", (int)nbosc / 2);
	gnuplot_set_xlabel(gp, "phase theta");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	gnuplot_plot_xy(gp, Nhisttheta, histthetafin, nbhist,"histogramme des valeurs finales des theta");	



	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 4 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 4 persist");
#endif
	if (nbK == 1) {

		/*Tracé de r, rmoyenRand et rayonmoyen en fonction du temps*/
//		gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//		gnuplot_cmd(gp, "set output 'rayon1.jpeg'");
		gnuplot_setstyle(gp, "dots");
		gnuplot_set_xlabel(gp, "t");
		gnuplot_set_ylabel(gp, "r");
		gnuplot_cmd(gp, "set yrange [-0.05:1.05]");
		sprintf(titre,"evolution de r(t) pour K = %f", K);
		gnuplot_plot_xy(gp, temps, rayon, nbsamples, titre);
		sprintf(titre,"evolution de rmoyen(t) pour K = %f", K);
		gnuplot_plot_xy(gp, temps, rayonmoyen, nbsamples, titre) ;
		sprintf(titre,"evolution de rmoy(t) pour K = %f", K);
		gnuplot_plot_xy(gp, temps, rayonmoyenRand, nbsamples, titre);



		/*Tracé de la phase des oscillateurs*/
		gp = gnuplot_init();
#if defined ( __APPLE__ )
		gnuplot_cmd(gp, "set terminal x11 5 persist");
#else
		gnuplot_cmd(gp, "set terminal wxt 5 persist");
#endif
//		gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//		gnuplot_cmd(gp, "set output 'angle.jpeg'");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_ylabel(gp, "phase theta");
		gnuplot_set_xlabel(gp, "numéro de l'oscillateur");
		gnuplot_plot_xy(gp, Nosc, omega, nbosc,"phase des oscillateurs");
	}


	else {

		/*Tracé de l'évolution de rstable et rinfini*/
//		gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//		gnuplot_cmd(gp, "set output 'rayon.jpeg'");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_xlabel(gp, "K");
		gnuplot_set_ylabel(gp, "r infini");
		gnuplot_cmd(gp, "set yrange [-0.05:1.05]");
		gnuplot_plot_xy(gp, Kvect, rayonstable, nbK, "evolution de rstable(K)");
		gnuplot_plot_xy(gp, Kvect, rayoninfini, nbK, "evolution de rinfini(K)");
	}

	printf("r = %f\npsi = %f\n",rayontemp,psitemp);



	/*Libération de la mémoire*/
	free(temps);
	free(omega);
	free(theta);
	free(rayon);
	free(psi);

	free(rayonmoyen);
	free(rayonmoyenRand);
	free(rayonstable);
	free(rayoninfini);
	free(Kvect);

	free(w);
	free(N);

	free(Tc);
	free(rayonCut);

	free(hist);
	free(histthetadeb);
	free(histthetafin);

	free(Nosc);
	free(Nhist);
	free(Nhisttheta);

	gnuplot_close(gp);
	return 0;

}

int meanField(double *theta, double *rayon, double *psi, int nbosc)
{
	int idxOsc;
	complex_d rComplex(0,0);
	for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
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
	return(omega + K * rayon * sin(psi - theta));
}

int logarithme(double Kc, double *Tc, double nbK, double Kmax, double *Kvect, double *rayonstable)
{
	int idxK;
	double ecartMax = 0;
	int idxprim = 0;
	int idxKc = 0;
	idxKc = (int)ceil(Kc * (nbK - 1) / Kmax) + 4; // le +4 est pour éviter d'être trop proche du seuil
	double *logk = (double *) malloc ((nbK-idxKc)*sizeof(double));
	double *logr = (double *) malloc ((nbK-idxKc)*sizeof(double));
	double *deltarstable = (double *) malloc ((nbK-idxKc)*sizeof(double));
	double *KvectCut = (double *) malloc ((nbK-idxKc)*sizeof(double));
	double *TcCut = (double *) malloc ((nbK-idxKc)*sizeof(double));

	for (idxK = idxKc ; idxK < nbK ; idxK++) {
		idxprim = idxK - idxKc;
		logk[idxprim] = log(1 - Kc / Kvect[idxK]);
		logr[idxprim] = log(rayonstable[idxK]);
		deltarstable[idxprim] = (sqrt(1 - Kc / Kvect[idxK]) - rayonstable[idxK]) / sqrt(1 - Kc / Kvect[idxK]);
		if (deltarstable[idxprim] > ecartMax)
			ecartMax = deltarstable[idxprim];
		KvectCut[idxprim] = Kvect[idxK];
		TcCut[idxprim] = Tc[idxK];
	}

	/*Régression linéaire*/
	double c0, c1, cov00, cov01, cov11, sumsq;
	size_t xstride = 1;
	size_t ystride = 1;
	size_t n = nbK - idxKc;
	double *fit = (double *) malloc (n*sizeof(double));
	gsl_fit_linear(logk, xstride, logr, ystride, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	for (idxK=0 ; idxK < nbK-idxKc ; idxK++) {
		fit[idxK] = c0 + c1 * logk[idxK];
	}



	/*Graphiques*/
	char titre[256];
	gnuplot_ctrl * gp;

	/*Tracé en log-log de rstable en fonction de (K-Kc)/K */
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 5 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 5 persist");
#endif
//	gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//	gnuplot_cmd(gp, "set output 'loglog.jpeg'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "log (K-Kc)/K");
	gnuplot_set_ylabel(gp, "log rstable");
//	gnuplot_cmd(gp, "set yrange [-0.05:10.05]");
//	gnuplot_plot_xy(gp, Kvect, invrayonstable, nbK, "evolution de rstable en fonction de K");
	sprintf(titre,"évolution de rstable en fonction de (K-Kc)/K en log-log");
	gnuplot_plot_xy(gp, logk, logr, (nbK-idxKc), titre);
	sprintf(titre,"régression linéaire");
	gnuplot_plot_xy(gp, logk, fit, (nbK-idxKc), titre);

	printf("écart entre simulation et théorie = %f\n",ecartMax);
	printf("pente du fit = %f\n", c1);
	printf("ordonnée à l'origine = %f\n", c0);



	/*Tracé de l'évolution du temps caractéractique en fonction de K*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 6 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 6 persist");
#endif
//	gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//	gnuplot_cmd(gp, "set output 'Tc.jpeg'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "K");
	gnuplot_set_ylabel(gp, "temps caractéristique");
	gnuplot_plot_xy(gp, KvectCut, TcCut, (nbK-idxKc), "évolution du temps caracteristique");



	/*Tracé de l'évolution de l'écart des r simulés et théoriques*/
	gp = gnuplot_init();
#if defined ( __APPLE__ )
	gnuplot_cmd(gp, "set terminal x11 7 persist");
#else
	gnuplot_cmd(gp, "set terminal wxt 7 persist");
#endif
//	gnuplot_cmd(gp, "set terminal jpeg enhanced color");
//	gnuplot_cmd(gp, "set output 'ecart.jpeg'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "K");
	gnuplot_set_ylabel(gp, "r théorique - r simulé");
	gnuplot_plot_xy(gp, KvectCut, deltarstable, (nbK-idxKc), "écart entre simulation et théorie");

	return 0;
}


