#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_histogram.h>

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
	size_t idxh;

	double rayontemp = 0, psitemp = 0;

	double *omega = (double *) malloc (nbosc*sizeof(double));
	double *theta = (double *) malloc (nbosc*sizeof(double));
	double *thetapoint = (double *) malloc (nbosc*sizeof(double));
	double K = 0;
	double Kmax = 1;
	double OMEGA = 0;
	double sigma = 0.1;
	double subcrit = 0;
	double *rayon = (double *) malloc (nbsamples*sizeof(double));
	double *psi = (double *) malloc (nbsamples*sizeof(double));

	double *rayonmoyenRand = (double *) malloc (nbsamples*sizeof(double)); /*moyenne sur les réalisations*/
//	double *invrayonstable = (double *) malloc (nbK*sizeof(double));
	double *Kvect = (double *) malloc (nbK*sizeof(double));

	/*Détermination de la distribution des pulsations propres*/
	int nbw = 1001;
	int idxw;
	double *w = (double *) malloc (nbw*sizeof(double));
	double *N = (double *) malloc (nbw*sizeof(double));

	/*Détermination de Kc et de alpha*/
	double Kc, alpha;
	Kc = 2 / (M_PI * gsl_ran_cauchy_pdf(subcrit, sigma));

	printf("Kc = %f\n", Kc);
	alpha = pow(Kc, 3) * (1 - 4 * pow(subcrit, 2) / (pow(sigma, 2) * (1 + pow(subcrit / sigma, 2))))/ (8 * pow(sigma, 3) * pow(1 + pow(subcrit / sigma, 2), 2));
	printf("alpha = %f\n", alpha);

	/*Déclaration des variables nécessaires à la détermination du rayon asymptotique rayonInfini*/
	int idxC = 0;
	double rayonMax = 0;
	double *Tc = (double *) malloc (nbK*sizeof(double));
	double *rayonCut = (double *) malloc (nbsamples*sizeof(double));
	double *rayoninfini = (double *) malloc (nbK*sizeof(double));
	double *rayonmoyen = (double *) malloc (nbsamples*sizeof(double));

	/*Déclaration des variables nécessaires pour Runge-Kutta 4*/

	double k1, k2, k3, k4;


	/*Définition des paramètres pour les histogrammes*/
	double wmax = 0.8;
	double thetamax = M_PI;
	size_t nhist = 100;
	int nbhist = (int) nhist;

	double *histw = (double *) malloc (nbhist*sizeof(double));
	double *histthetadeb = (double *) malloc (nbhist*sizeof(double));
	double *histthetafin = (double *) malloc (nbhist*sizeof(double));
	double *histthetapointdeb = (double *) malloc (nbhist*sizeof(double));
	double *histthetapointfin = (double *) malloc (nbhist*sizeof(double));

	double thetapointmin = -0.1;
	double thetapointmax = 0.1;

	gsl_histogram *htheta = gsl_histogram_alloc(nhist);
	gsl_histogram *hthetapoint = gsl_histogram_alloc(nhist);
	gsl_histogram_set_ranges_uniform (htheta, 0, 2 * thetamax);
	gsl_histogram_set_ranges_uniform (hthetapoint, thetapointmin, thetapointmax);

	FILE* fichier = NULL;

	/*Détermination des pulsations propres des oscillateurs*/
	gsl_rng * r;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	gsl_rng_default_seed = rand();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
//		printf ("seed = %lu\n", gsl_rng_default_seed);

	for (idxOsc = 0 ; idxOsc < (nbosc / 2) ; idxOsc++) {
//			omega[idxOsc] = OMEGA + gsl_ran_gaussian(r,sigma) + subcrit;
		omega[idxOsc] = OMEGA + gsl_ran_cauchy(r,sigma) + subcrit;
	}

	for (idxOsc = (nbosc / 2) ; idxOsc < nbosc ; idxOsc++) {
//			omega[idxOsc] = OMEGA + gsl_ran_gaussian(r,sigma) - subcrit;
		omega[idxOsc] = OMEGA + gsl_ran_cauchy(r,sigma) - subcrit;
	}


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

		/*Initialisation de rayonmoyenRand*/
		for (idxTime=0 ; idxTime < nbsamples ; idxTime++) {
			rayonmoyenRand[idxTime] = 0;
		}


		/*Boucle sur les réalisations*/
		for (idxRand = 0 ; idxRand < nbrand ; idxRand++) {


			/*Histogramme des pulsations*/
			if (idxRand == 0) {
				gsl_histogram *h = gsl_histogram_alloc(nhist);
				gsl_histogram_set_ranges_uniform (h, -wmax, wmax);

				for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
					gsl_histogram_increment (h, omega[idxOsc]);
				}

				for (idxh = 0 ; idxh < nbhist ; idxh++) {
					histw[idxh] = gsl_histogram_get (h, idxh);
				}
				gsl_histogram_free (h);
			}

			/*Initialisation de theta, de rayon et de psi*/
			for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
				theta[idxOsc] = modulo(rand()/10000.0, 2 * M_PI);
			}

			for (idxTime = 0 ; idxTime < nbsamples ; idxTime++) {
				rayon[idxTime] = 0;
				psi[idxTime] = 0;
			}

			meanField(theta , &rayontemp, &psitemp, nbosc);
			rayon[0] = rayontemp;
			psi[0] = psitemp;

			/*Corps du programme*/
			/*Histogramme des valeurs initiales des theta*/
			if (idxK == nbK - 1 && idxRand == 0) {
				for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
					theta[idxOsc] = modulo(theta[idxOsc], 2 * M_PI);
					gsl_histogram_increment (htheta, theta[idxOsc]);
				}

				for (idxh = 0 ; idxh < nbhist ; idxh++) {
					histthetadeb[idxh] = gsl_histogram_get (htheta, idxh);
				}
			}
			for (idxTime = 1 ; idxTime < nbsamples ; idxTime++) {
				temps[idxTime] = idxTime * deltaT;

				for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
					thetapoint[idxOsc]=theta[idxOsc];
					/*Détermination des theta par la méthode de runge-kutta d'ordre 4*/
					k1 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc]);
					k2 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k1 / 2.0);
					k3 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k2 / 2.0);
					k4 = deltaT * kuramoto(omega[idxOsc], K, psi[idxTime - 1], rayon[idxTime - 1], theta[idxOsc] + deltaT * k3);
					theta[idxOsc] = theta[idxOsc] + (k1 + 2 * k2 + 2 * k3 + k4) * deltaT / 6.0;
					thetapoint[idxOsc] = (theta[idxOsc] - thetapoint[idxOsc]) / deltaT;
				}

				/*Détermination de r et psi par l'approche champ moyen*/
				meanField(theta , &rayontemp, &psitemp, nbosc);
				rayon[idxTime] = rayontemp;
				psi[idxTime] = psitemp;

				/*Détermination de rayonmoyenRand*/
				rayonmoyenRand[idxTime] += rayontemp / nbrand;
	
				
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == nbsamples - 1) {
					fichier = fopen("thetapointfin.csv", "w");
				
				    if (fichier != NULL) {
						for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {	
				   			fprintf(fichier, "%d;%f\n", idxOsc, thetapoint[idxOsc]);
						}
				        fclose(fichier);
				    }
					else {
				        printf("Impossible d'ecrire dans le fichier");
				    }
				}
				

				/*Histogramme des valeurs finales des theta*/
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == nbsamples - 1) {
					gsl_histogram_reset (htheta);
					for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
						theta[idxOsc] = modulo(theta[idxOsc], 2 * M_PI);
						gsl_histogram_increment (htheta, theta[idxOsc]);
					}

					for (idxh = 0 ; idxh < nbhist ; idxh++) {
						histthetafin[idxh] = gsl_histogram_get (htheta, idxh);
					}
					gsl_histogram_free (htheta);
				}

				/*Histogramme des valeurs initiales des thetapoint*/
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == 1) {
					for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
						gsl_histogram_increment (hthetapoint, thetapoint[idxOsc]);
					}

					for (idxh = 0 ; idxh < nbhist ; idxh++) {
						histthetapointdeb[idxh] = gsl_histogram_get (hthetapoint, idxh);
					}
				}

				/*Histogramme des valeurs finales des thetapoint*/
				if (idxK == nbK - 1 && idxRand == 0 && idxTime == nbsamples - 1) {
					gsl_histogram_reset (hthetapoint);
					for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
						gsl_histogram_increment (hthetapoint, thetapoint[idxOsc]);
					}

					for (idxh = 0 ; idxh < nbhist ; idxh++) {
						histthetapointfin[idxh] = gsl_histogram_get (hthetapoint, idxh);
					}
					gsl_histogram_free (hthetapoint);
				}


			}

			/*Définition de rayonMax, le maximum de r*/
			for (idxTime = 0 ; idxTime < nbsamples ; idxTime++) {
				if (rayonmoyenRand[idxTime] > rayonMax)
					rayonMax = rayonmoyenRand[idxTime];
			}

			/*Définition du temps critique comme étant le temps pour lequel on atteint 90% de la valeur maximale*/
			idxTime = 0;
			idxC = 0;
			while (rayonmoyenRand[idxTime] < 0.90 * rayonMax)
				idxTime++;

			idxC = idxTime;
			if (idxC < nbsamples)
				Tc[idxK] += idxC*deltaT / nbrand;
		}

		/*Fin de la boucle sur les réalisations*/

//		invrayonstable[idxK] = 1 / (1 - rayonmoyen[idxK] * rayonmoyen[idxK]);	//pour formule 4.7


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


		/*Détermination de rayonmoyen, comme moyenne temporelle glissante sur les nbmoyglissante derniers temps*/
		if (idxK == nbK - 1) {
			int nbmoyglissante = 100;
			rayonmoyen[0] = rayon[0];
			for (idxTime = 1 ; idxTime < nbmoyglissante ; idxTime++) {
				rayonmoyen[idxTime] = rayonmoyen[idxTime - 1] + rayon[idxTime];
			}
			for (idxTime = nbmoyglissante ; idxTime < nbsamples ; idxTime++) {
				rayonmoyen[idxTime] = rayonmoyen[idxTime - 1] + rayon[idxTime] - rayon[idxTime - nbmoyglissante];
			}
			for (idxTime = 0 ; idxTime < nbsamples ; idxTime++) {
				rayonmoyen[idxTime] /= nbmoyglissante;
			}
		}
	}


	/*Définition pour les logarithmes*/
	if (nbK > 4 && Kc < Kvect[nbK - 1]) {
		logarithme(Kc, Tc, nbK, Kmax, Kvect, rayoninfini);
	}

	/*Détermination de la distribution des pulsations propres*/
	for (idxw=0 ; idxw < nbw ; idxw++) {
		w[idxw] = (-(nbw - 1) / 2 + idxw) * wmax * 2 / nbw;
		N[idxw] = gsl_ran_cauchy_pdf(w[idxw] - subcrit, sigma) / 2 + gsl_ran_cauchy_pdf(w[idxw] + subcrit, sigma) / 2;
	}

	/*Détermination des vecteurs des abscisses pour les distributions et les histogrammes*/
	double *Nosc = (double *) malloc (nbosc*sizeof(double));
	double *Nhist = (double *) malloc (nbhist*sizeof(double));
	double *Nhistthetapoint = (double *) malloc (nbhist*sizeof(double));
	double *Nhisttheta = (double *) malloc (nbhist*sizeof(double));

	for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {
		Nosc[idxOsc] = idxOsc;
	}
	for (idxOsc = 0 ; idxOsc < nbhist ; idxOsc++) {
		Nhist[idxOsc] = idxOsc * 2.0 * wmax / nbhist - wmax;
		Nhisttheta[idxOsc] = idxOsc * 2.0 * thetamax / nbhist; //XXX
		Nhistthetapoint[idxOsc] = idxOsc * (thetapointmax - thetapointmin) / nbhist - thetapointmax;
	}




	/*Définitions pour le graphique*/
	gnuplot_ctrl * gp;
	char titre[256];


	/*Tracé de la distribution théorique des pulsations*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'distribution.ps'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "pulsation w");
	gnuplot_set_ylabel(gp, "distribution g(w)");
	gnuplot_plot_xy(gp, w, N, nbw,"distribution des pulsations propres des oscillateurs");


	/*Tracé de l'histogramme des valeurs des pulsations*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'histo-pulsations.ps'");
	gnuplot_setstyle(gp, "steps");
	gnuplot_set_xlabel(gp, "pulsation w");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	gnuplot_plot_xy(gp, Nhist, histw, nbhist,"histogramme des valeurs des pulsations");



	/*Tracé de l'histogramme des valeurs initiales des theta*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'histo-theta-init.ps'");
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set xrange [-0.2:%f]", 2 * thetamax + 0.2);
	gnuplot_cmd(gp, "set yrange [0:100]");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	gnuplot_set_xlabel(gp, "phase theta");
	sprintf(titre,"histogramme des valeurs initiales des theta pour K = %f", K);
	gnuplot_plot_xy(gp, Nhisttheta, histthetadeb, nbhist, titre);


	/*Tracé de l'histogramme des valeurs finales des theta*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'histo-theta-fin.ps'");
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set xrange [-0.2:%f]", 2 * thetamax + 0.2);
	gnuplot_cmd(gp, "set yrange [0:100]");
	gnuplot_set_xlabel(gp, "phase theta");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	sprintf(titre,"histogramme des valeurs finales des theta pour K = %f", K);
	gnuplot_plot_xy(gp, Nhisttheta, histthetafin, nbhist, titre);

	/*Tracé de l'histogramme des valeurs initiales des thetapoint*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'histo-thetapoint-deb.ps'");
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set xrange [%f:%f]", thetapointmin, thetapointmax);
	gnuplot_cmd(gp, "set yrange [0:%d]", nbosc);
	gnuplot_set_xlabel(gp, "phase theta");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	sprintf(titre,"histogramme des valeurs initiales des thetapoint pour K = %f", K);
	gnuplot_plot_xy(gp, Nhistthetapoint, histthetapointdeb, nbhist, titre);

	/*Tracé de l'histogramme des valeurs finales des thetapoint*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'histo-thetapoint-fin.ps'");
	gnuplot_setstyle(gp, "steps");
	gnuplot_cmd(gp, "set xrange [%f:%f]", thetapointmin, thetapointmax);
	gnuplot_cmd(gp, "set yrange [0:%d]", nbosc);
	gnuplot_set_xlabel(gp, "phase theta");
	gnuplot_set_ylabel(gp, "Nombre d'oscillateurs");
	sprintf(titre,"histogramme des valeurs finales des thetapoint pour K = %f", K);
	gnuplot_plot_xy(gp, Nhistthetapoint, histthetapointfin, nbhist, titre);


	/*fichier = fopen("histthetadeb.csv", "w");

    if (fichier != NULL) {
		for (idxh = 0 ; idxh < nbhist ; idxh++) {	
   			fprintf(fichier, "%f;%f\n", Nhisttheta[idxh], histthetadeb[idxh]);
		}
        fclose(fichier);
    }
	else {
        printf("Impossible d'ecrire dans le fichier");
    }
	
	fichier = fopen("histthetafin.csv", "w");

    if (fichier != NULL) {
		for (idxh = 0 ; idxh < nbhist ; idxh++) {	
   			fprintf(fichier, "%f;%f\n", Nhisttheta[idxh], histthetafin[idxh]);
		}
        fclose(fichier);
    }
	else {
        printf("Impossible d'ecrire dans le fichier");
    }

	fichier = fopen("thetafin.csv", "w");

    if (fichier != NULL) {
		for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++) {	
   			fprintf(fichier, "%d;%f\n", idxOsc, theta[idxOsc]);
		}
        fclose(fichier);
    }
	else {
        printf("Impossible d'ecrire dans le fichier");
    }
	*/


	
	if (nbK == 1) {

		/*Tracé de r, rmoyenRand et rayonmoyen en fonction du temps*/
		gp = gnuplot_init();
		gnuplot_cmd(gp, "set terminal postscript enhanced color");
		gnuplot_cmd(gp, "set output \"rayon1.ps\"");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_xlabel(gp, "t");
		gnuplot_set_ylabel(gp, "r");
		gnuplot_cmd(gp, "set yrange [-0.05:1.05]");
//		sprintf(titre,"evolution de r(t) pour K = %f", K);
//		gnuplot_plot_xy(gp, temps, rayon, nbsamples, titre);
//		sprintf(titre,"evolution de rayonmoyen(t) pour K = %f", K);
//		gnuplot_plot_xy(gp, temps, rayonmoyen, nbsamples, titre) ;
		sprintf(titre,"evolution de rmoy(t) pour K = %f", K);
		gnuplot_plot_xy(gp, temps, rayonmoyenRand, nbsamples, titre);

		/*Tracé de psi*/
		gp = gnuplot_init();
		gnuplot_cmd(gp, "set terminal postscript enhanced color");
		gnuplot_cmd(gp, "set output \"psi1.ps\"");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_xlabel(gp, "t");
		gnuplot_set_ylabel(gp, "psi");
		gnuplot_cmd(gp, "set yrange [-3.2:3.2]");
		sprintf(titre,"evolution de psi(t) pour K = %f", K);
		gnuplot_plot_xy(gp, temps, psi, nbsamples, titre);

		/*Tracé de la phase des oscillateurs*/
		gp = gnuplot_init();
		gnuplot_cmd(gp, "set terminal postscript enhanced color");
		gnuplot_cmd(gp, "set output \"angle.ps\"");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_ylabel(gp, "phase theta");
		gnuplot_set_xlabel(gp, "numero de l'oscillateur");
		gnuplot_plot_xy(gp, Nosc, omega, nbosc,"phase des oscillateurs");
	}


	else {

		/*Tracé de l'évolution de rstable et rinfini*/
		gp = gnuplot_init();
		gnuplot_cmd(gp, "set terminal postscript enhanced color");
		gnuplot_cmd(gp, "set output 'rayon.ps'");
		gnuplot_setstyle(gp, "lines");
		gnuplot_set_xlabel(gp, "K");
		gnuplot_set_ylabel(gp, "r infini");
		gnuplot_cmd(gp, "set yrange [-0.05:1.05]");
		gnuplot_plot_xy(gp, Kvect, rayoninfini, nbK, "evolution de rinfini(K)");
	}

	printf("r = %f\npsi = %f\n",rayontemp,psitemp);
	printf("K.r = %f\n", K * rayoninfini[nbK - 1]);


	/*Libération de la mémoire*/
	free(temps);
	free(omega);
	free(theta);
	free(rayon);
	free(psi);

	free(rayonmoyen);
	free(rayonmoyenRand);
	free(rayoninfini);
	free(Kvect);

	free(w);
	free(N);

	free(Tc);
	free(rayonCut);

	free(histw);
	free(histthetadeb);
	free(histthetafin);
	free(histthetapointdeb);
	free(histthetapointfin);

	free(Nosc);
	free(Nhist);
	free(Nhisttheta);
	free(Nhistthetapoint);

	gnuplot_close(gp);
	return 0;

}

int meanField(double *theta, double *rayon, double *psi, int nbosc)
{
	int idxOsc;
	gsl_complex rComplex;
	gsl_complex rComplexTemp;
	rComplex = gsl_complex_rect(0,0);
	for (idxOsc = 0 ; idxOsc < nbosc ; idxOsc++)
	{
		rComplexTemp = gsl_complex_polar(1,theta[idxOsc]);
		rComplex = gsl_complex_add(rComplex,rComplexTemp);
	}
	rComplex = gsl_complex_div_real(rComplex,nbosc);

	*rayon = gsl_complex_abs(rComplex);
	*psi = gsl_complex_arg(rComplex);

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
		TcCut[idxprim] = log(Tc[idxK]);
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
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'loglog.ps'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "k");
	gnuplot_set_ylabel(gp, "rayon infini");
//	gnuplot_cmd(gp, "set yrange [-0.05:10.05]");
//	gnuplot_plot_xy(gp, Kvect, invrayonstable, nbK, "evolution de rstable en fonction de K");
	sprintf(titre,"evolution du rayon infini");
	gnuplot_plot_xy(gp, logk, logr, (nbK-idxKc), titre);
	sprintf(titre,"regression lineaire");
	gnuplot_plot_xy(gp, logk, fit, (nbK-idxKc), titre);

	printf("écart entre simulation et théorie = %f\n",ecartMax);
	printf("pente du fit = %f\n", c1);
	printf("coefficient de corrélation = %f\n", gsl_stats_correlation(logk, xstride, logr, ystride, n));
	printf("alpha simulé = %f\n", pow(10, -2 * c0));



	/*Tracé de l'évolution du temps caractéractique en fonction de K*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'Tc.ps'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "K");
	gnuplot_set_ylabel(gp, "temps caracteristique");
	gnuplot_plot_xy(gp, KvectCut, TcCut, (nbK-idxKc), "Evolution du temps caracteristique");



	/*Tracé de l'évolution de l'écart des r simulés et théoriques*/
	gp = gnuplot_init();
	gnuplot_cmd(gp, "set terminal postscript enhanced color");
	gnuplot_cmd(gp, "set output 'ecart.ps'");
	gnuplot_setstyle(gp, "lines");
	gnuplot_set_xlabel(gp, "K");
	gnuplot_set_ylabel(gp, "r théorique - r simulé");
	gnuplot_plot_xy(gp, KvectCut, deltarstable, (nbK-idxKc), "Ecart entre simulation et théorie");

	return 0;
}

double modulo(double num,double deno)
{
	deno = abs(deno);
	while (num > 0) {
		num -= deno;
	}
	while (num < 0) {
		num += deno;
	}
	return num;
}
