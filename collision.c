#include <stdio.h>
#include <math.h>
#include "collision.h"
#include "config.h"
#include "utils.h"

double collision_position_time(double *rf, double *rb,
                               double *vs, double *vb,
                               int *shift,
                               double *rfc, double *rbc,
                               double *nc, struct Config *cf) {
	//printf("collision_sphere_exact\n");
	// from lammps
	double vs_dot_vs,vb_dot_vb,vs_dot_vb;
	double vs_dot_rb,vb_dot_rf,vs_dot_rf,vb_dot_rb;
	double rf_dot_rf,rb_dot_rb,rf_dot_rb;
	double a,b,c,scale;
	double rbp[3];
	double dtt;

	rbp[0] = rb[0] + shift[0]*cf->lx;
	rbp[1] = rb[1] + shift[1]*cf->ly;
	rbp[2] = rb[2] + shift[2]*cf->lz;

	vs_dot_vs = vs[0]*vs[0] + vs[1]*vs[1] + vs[2]*vs[2];
	vb_dot_vb = vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2];
	vs_dot_vb = vs[0]*vb[0] + vs[1]*vb[1] + vs[2]*vb[2];

	vs_dot_rb = vs[0]*rbp[0] + vs[1]*rbp[1] + vs[2]*rbp[2];
	vb_dot_rf = vb[0]*rf[0]  + vb[1]*rf[1]  + vb[2]*rf[2];
	vs_dot_rf = vs[0]*rf[0]  + vs[1]*rf[1]  + vs[2]*rf[2];
	vb_dot_rb = vb[0]*rbp[0] + vb[1]*rbp[1] + vb[2]*rbp[2];

	rf_dot_rf = rf[0]*rf[0]   + rf[1]*rf[1]   + rf[2]*rf[2];
	rb_dot_rb = rbp[0]*rbp[0] + rbp[1]*rbp[1] + rbp[2]*rbp[2];
	rf_dot_rb = rf[0]*rbp[0]  + rf[1]*rbp[1]  + rf[2]*rbp[2];

	a = vs_dot_vs + vb_dot_vb - 2.0*vs_dot_vb;
	b = 2.0 * (vs_dot_rb + vb_dot_rf - vs_dot_rf - vb_dot_rb);
	c = rf_dot_rf + rb_dot_rb - 2.0*rf_dot_rb - cf->rb*cf->rb;

	dtt = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);

	rfc[0] = rf[0] - dtt*vs[0];
	rfc[1] = rf[1] - dtt*vs[1];
	rfc[2] = rf[2] - dtt*vs[2];

	rbc[0] = rbp[0] - dtt*vb[0];
	rbc[1] = rbp[1] - dtt*vb[1];
	rbc[2] = rbp[2] - dtt*vb[2];

	nc[0] = rfc[0] - rbc[0];
	nc[1] = rfc[1] - rbc[1];
	nc[2] = rfc[2] - rbc[2];
	scale = 1.0/sqrt(nc[0]*nc[0] + nc[1]*nc[1] + nc[2]*nc[2]);
	nc[0] *= scale;
	nc[1] *= scale;
	nc[2] *= scale;

	return dtt;
}

void collision_velocities(
	double *rfc, double *vf,
	double *rbc, double *vb,
	double *wb, double *ob,
	double *nc,
	double *vfn,
	int it,
	struct Config *cf) {
	
	//printf("noslip\n");
	
	double vf_dot_nc;
	double ob_dot_nc;
	double scale,r1,r2,vnmag,vtmag1,vtmag2;
	double tv1[3]; // tangent to big particle surface, in direction of the fluid particle velocity
	double tv2[3];
	double to1[3]; // tangent to big particle surface, in direction formed between normal and orientation 
	double to2[3];
	double lamda, sigma, dmax, vmax, vmaxsq;
	double vth;
	
	lamda = cf->dt*sqrt(cf->Tf/cf->mf);
	sigma = lamda/cf->dt;
	dmax = 4.0*lamda;
	vmax = dmax/cf->dt;
  	vmaxsq = vmax*vmax;

	vf_dot_nc = vf[0]*nc[0] + vf[1]*nc[1] + vf[2]*nc[2];
	tv1[0] = vf[0] - vf_dot_nc*nc[0];
	tv1[1] = vf[1] - vf_dot_nc*nc[1];
	tv1[2] = vf[2] - vf_dot_nc*nc[2];
  	scale = 1.0/sqrt(tv1[0]*tv1[0] + tv1[1]*tv1[1] + tv1[2]*tv1[2]);
  	tv1[0] *= scale;
  	tv1[1] *= scale;
  	tv1[2] *= scale;

  	tv2[0] = nc[1]*tv1[2] - nc[2]*tv1[1];
  	tv2[1] = nc[2]*tv1[0] - nc[0]*tv1[2];
  	tv2[2] = nc[0]*tv1[1] - nc[1]*tv1[0];

  	ob_dot_nc = ob[0]*nc[0] + ob[1]*nc[1] + ob[2]*nc[2];
  	to1[0] = ob[0] - (ob_dot_nc)*nc[0];
  	to1[1] = ob[1] - (ob_dot_nc)*nc[1];
  	to1[2] = ob[2] - (ob_dot_nc)*nc[2];
  	scale = 1.0/sqrt(to1[0]*to1[0] + to1[1]*to1[1] + to1[2]*to1[2]);
  	to1[0] *= scale;
  	to1[1] *= scale;
  	to1[2] *= scale;

	to2[0] = nc[1]*to1[2] - nc[2]*to1[1];
  	to2[1] = nc[2]*to1[0] - nc[0]*to1[2];
  	to2[2] = nc[0]*to1[1] - nc[1]*to1[0];

  	while (1) {
    	r1 = sigma * gaussian();
    	r2 = sigma * gaussian();
    	vnmag = sqrt(r1*r1 + r2*r2);
    	vtmag1 = sigma * gaussian();
    	vtmag2 = sigma * gaussian();
    	if (vnmag*vnmag + vtmag1*vtmag1 + vtmag2*vtmag2 <= vmaxsq) break;
  	}

	double cb1 = 0.0;
	double cb2 = 0.0;

	double dx = rfc[0]-rbc[0];
	double dy = rfc[1]-rbc[1];
	double dz = rfc[2]-rbc[2];

	double theta;
  	double geom = (dx*ob[0] + dy*ob[1] + dz*ob[2])/cf->rb;
  	if (geom > 1.0) {
  		printf("Collision angle overflow (+)\n");
  		printf("dx=%f dy=%f dz=%f\n",dx,dy,dz);
  		printf("ob=(%.16f %.16f %.16f)\n",ob[0],ob[1],ob[2]);
  		printf("geom=%.16f\n",geom);
  		geom = 1.0;
  	} else if (geom < -1.0) {
  		printf("Collision angle overflow (-)\n");
  		printf("dx=%f dy=%f dz=%f\n",dx,dy,dz);
  		printf("ob=(%.16f %.16f %.16f)\n",ob[0],ob[1],ob[2]);
  		printf("geom=%.16f\n",geom);
  		geom = -1.0;
  	}
  	theta = acos(geom);

	vth = (3.0/3.0)*(cb1 + cb2*cos(theta))*sin(theta);

	// velocity given by no slip boundary condition

  	vfn[0] = vnmag*nc[0] + (vtmag1)*tv1[0] + vtmag2*tv2[0];
  	vfn[1] = vnmag*nc[1] + (vtmag1)*tv1[1] + vtmag2*tv2[1];
  	vfn[2] = vnmag*nc[2] + (vtmag1)*tv1[2] + vtmag2*tv2[2];

  	// velocity added by squirmer

  	vfn[0] += (vth)*to1[0];
  	vfn[1] += (vth)*to1[1];
  	vfn[2] += (vth)*to1[2];

  	// add velocity of big particle

	vfn[0] += vb[0] + wb[1]*dz - wb[2]*dy;
	vfn[1] += vb[1] + wb[2]*dx - wb[0]*dz;
	vfn[2] += vb[2] + wb[0]*dy - wb[1]*dx;

	
	FILE * f;
	f = fopen("collisions.txt","a");
	fprintf(f,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
		it,
		rbc[0],rbc[1],rbc[2],
		rfc[0],rfc[1],rfc[2],
		vb[0],vb[1],vb[2],
		vf[0],vf[1],vf[2],
		vfn[0],vfn[1],vfn[2],
		ob[0],ob[1],ob[2],
		theta,vth,
		nc[0],nc[1],nc[2]);
	fclose(f);
	
}

void collision_forces(double *vfold, double *vfnew,
					  double *rfc, double *rbc,
					  double *fb, double *tb, struct Config *cf) {

	//printf("force_torque\n");
	double dpdt[3],xf_xb[3];
	double factor = cf->mf / cf->dt;

	dpdt[0] = factor * (vfnew[0] - vfold[0]);
	dpdt[1] = factor * (vfnew[1] - vfold[1]);
	dpdt[2] = factor * (vfnew[2] - vfold[2]);

	fb[0] -= dpdt[0];
	fb[1] -= dpdt[1];
	fb[2] -= dpdt[2];

	xf_xb[0] = rfc[0]-rbc[0];
	xf_xb[1] = rfc[1]-rbc[1];
	xf_xb[2] = rfc[2]-rbc[2];

	// t = r x f, this is the torque with respect to center of particle
	tb[0] -= xf_xb[1]*dpdt[2] - xf_xb[2]*dpdt[1];
	tb[1] -= xf_xb[2]*dpdt[0] - xf_xb[0]*dpdt[2];
	tb[2] -= xf_xb[0]*dpdt[1] - xf_xb[1]*dpdt[0];
}