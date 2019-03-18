#pragma once

struct FParticle {
	double r[3];	// position
	double v[3];	// velocity
};

struct BParticle {
	double r[3];	// position
	double v[3];	// velocity
	double a[3];	// acceleration
	double f[3];	// force
	double o[3];	// orientation
	double aa[3];	// angular acceleration
	double w[3];	// angular velocity
	double t[3];	// torque
};