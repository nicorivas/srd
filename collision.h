#pragma once
#include "config.h"

double collision_position_time(double *rs, double *rb,
                               double *vs, double *vb,
                               int *shift,
                               double *rsc, double *rbc,
                               double *norm, struct Config *cf);

void collision_velocities(double *vf, double *vb,
						  double *rb, double *wb, double *ob,
						  double *rc, double *nc,
						  double *vfn, int it, struct Config *cf);

void collision_forces(double *vsold, double *vsnew,
					  double *xs, double *xb,
					  double *fb, double *tb, struct Config *cf);