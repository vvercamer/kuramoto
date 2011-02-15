#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "gnuplot.h"

static int
lserial = 0;

GNUplot::GNUplot ()

{
	char cmd[256];

	lserial++;
	serial = lserial;

	snprintf(cmd, sizeof(cmd), "gnuplot -persist");
	gnuplotpipe = popen(cmd, "w");
	if (!gnuplotpipe) {
		throw("gnuplot not found !");
	}

	uid = (int)getuid();
} 

GNUplot::~GNUplot () 

{
	fprintf(gnuplotpipe, "exit\n");
	pclose(gnuplotpipe);
} 

void 
GNUplot::operator () (const char * command) 

{
	fprintf(gnuplotpipe, "%s\n", command);
	fflush(gnuplotpipe);
}

void
GNUplot::draw (	double * x,
		double * y1, 
		double * y2,
		int n)

{
	FILE * fp;
	int i;
	double v, max = 0.0, tmin = 0.0, tmax = 0.0;
	char cmd[256];

	for (i = 0; i < n; i++) {
		v = fabs(y1[i]);
		if (v > max)
			max = v;
		v = fabs(y2[i]);
		if (v > max)
			max = v;
		if (i == 0 || x[i] > tmax)
			tmax = x[i];
		if (i == 0 || x[i] < tmin)
			tmin = x[i];
	}

	snprintf(cmd, sizeof(cmd), "/tmp/%d-1-%d.dat", uid, serial); 
	fp = fopen(cmd, "w");
	if (fp != NULL) {
		for (i = 0; i < n; i++) 
			fprintf(fp, "%g %g\n", x[i], y1[i]);
		fclose(fp);
	}

	snprintf(cmd, sizeof(cmd), "/tmp/%d-2-%d.dat", uid, serial); 
	fp = fopen(cmd, "w");
	if (fp != NULL) {
		for (i = 0; i < n; i++) 
			fprintf(fp, "%g %g\n", x[i], y2[i]);
		fclose(fp);
	}

	snprintf(cmd, sizeof(cmd), 
		"plot [%g:%g] [%g:%g] '/tmp/%d-1-%d.dat' w l, '/tmp/%d-2-%d.dat' w l",
		tmin, tmax, -max, max, uid, serial, uid, serial);
	(*this)(cmd);
}

void
GNUplot::draw (	double * x,
		double * y, 
		int n)

{
	FILE * fp;
	int i;
	double v, min = 0.0, max = 0.0, tmin = 0.0, tmax = 0.0;
	char cmd[256];

	min = max = y[0];
	for (i = 0; i < n; i++) {
		v = y[i];
		if (v > max)
			max = v;
		if (v < min)
			min = v;
		if (i == 0 || x[i] > tmax)
			tmax = x[i];
		if (i == 0 || x[i] < tmin)
			tmin = x[i];
	}

	snprintf(cmd, sizeof(cmd), "/tmp/%d-1-%d.dat", uid, serial); 
	fp = fopen(cmd, "w");
	if (fp != NULL) {
		for (i = 0; i < n; i++) 
			fprintf(fp, "%g %g\n", x[i], y[i]);
		fclose(fp);
	}

	snprintf(cmd, sizeof(cmd), 
		"plot [%g:%g] [%g:%g] '/tmp/%d-1-%d.dat' w l",
		tmin, tmax, min, max, uid, serial);
	(*this)(cmd);
}

