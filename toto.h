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
int meanField(double *theta, double *rayon, double *psi, int nboscs);

double kuramoto(double omega, double K, double psi, double rayon, double theta);

double gaussianRand(void);
#endif
