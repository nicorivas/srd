#pragma once

#include "globals.h"

double randf();
double gaussian();
int mod(int a, int b);
double nmodf(double a, double b);
void random_unitary_vector(float *v);
void matrix_inverse(double m[DIM][DIM], double minv[DIM][DIM]);