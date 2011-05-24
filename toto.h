#ifndef _TOTO_H_INCLUDED
#define _TOTO_H_INCLUDED
/*
typedef struct KuramotoStruct
{
	double *omega;
	double K;
	double *rayon;
	double *psi;
} KuramotoS;
*/
int meanField(double *theta, double *rayon, double *psi, int nbsamples);

double kuramoto(double omega, double K, double psi, double rayon, double theta);

double gaussianRand(void);

int logarithme(double Kc, double *Tc, double nbK, double Kmax, double *Kvect, double *rayonstable);

double modulo(double num, double deno);
#endif
