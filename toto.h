#ifndef _TOTO_H_INCLUDED
#define _TOTO_H_INCLUDED

typedef struct KuramotoStruct KuramotoStruct;
struct KuramotoStruct
{
	double *omega;
	double K;
	double *rayon;
	double *psi;
};

int meanField(double *theta, double *rayon, double *psi, int nboscs);

#endif
