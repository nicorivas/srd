#include <stdlib.h>
#include <math.h>
#include "utils.h"

double randf() {
	return (double)rand()/(double)(RAND_MAX);
}

double gaussian() {
	double first,v1,v2,rsq,fac;
	do {
		v1 = 2.0*randf()-1.0;
		v2 = 2.0*randf()-1.0;
		rsq = v1*v1 + v2*v2;
	} while ((rsq >= 1.0) || (rsq == 0.0));
	return v1*sqrt(-2.0*log(rsq)/rsq);
}

int mod(int a, int b) {
	int ret = a % b;
	if (ret < 0)
		ret+=b;
	return ret;
}

double nmodf(double a, double b) {
	double ret = modf(a,&b);
	if (ret < 0)
		ret += b;
	return ret;
}

void random_unitary_vector(float *v) {
	double r1 = randf();
	double r2 = randf();
	double phi = 2.0*M_PI*r1;
	double rho = 2.0*r2-1.0;
	v[0] = cos(phi)*sqrt(1.0-rho*rho);
	v[1] = sin(phi)*sqrt(1.0-rho*rho);
	v[2] = rho;
}

void matrix_inverse(double m[DIM][DIM], double minv[DIM][DIM]) {
	// computes the inverse of a matrix m
	double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
	             m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
	             m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

	double invdet = 1.0 / det;

	minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
	minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
	minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
	minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
	minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
	minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
	minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
	minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
	minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}