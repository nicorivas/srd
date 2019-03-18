#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "config.h"
#include "collision.h"
#include "utils.h"
#include "output.h"
#include "particles.h"
#include "cells.h"

#define DIM 3
//#define AMC
//#define DEBUG_LOOP

double sqrttfmf;

double distance(double *p1, double *p2, int *shift, struct Config *cf) {
	int dx, dy, dz, dmin, d;
	dmin = cf->lx;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			for (dz=-1; dz<=1; dz++) {
				d = (p1[0] - p2[0] - cf->lx*dx)*(p1[0] - p2[0] - cf->lx*dx) +
	 		        (p1[1] - p2[1] - cf->ly*dy)*(p1[1] - p2[1] - cf->ly*dy) +
	 		        (p1[2] - p2[2] - cf->lz*dz)*(p1[2] - p2[2] - cf->lz*dz);
				if (d < dmin) {
					dmin = d;
					shift[0] = dx;
					shift[1] = dy;
					shift[2] = dz;
				}
			}
		}
	}
	return sqrt(dmin);
}

void update_srd(double dtt, double *xscoll, double *vsnew,
               double *xs, double *vs, struct Config *cf) {
	//printf("update_srd\n");
	vs[0] = vsnew[0];
	vs[1] = vsnew[1];
	vs[2] = vsnew[2];

	xs[0] = xscoll[0] + dtt*vsnew[0];
	xs[1] = xscoll[1] + dtt*vsnew[1];
	xs[2] = xscoll[2] + dtt*vsnew[2];

	if (xs[0] < 0.0)
		xs[0] += cf->lx;
	else if (xs[0] >= cf->lx)
		xs[0] -= cf->lx;

	if (xs[1] < 0.0)
		xs[1] += cf->ly;
	else if (xs[1] >= cf->ly)
		xs[1] -= cf->ly;

	if (xs[2] < 0.0)
		xs[2] += cf->lz;
	else if (xs[2] >= cf->lz)
		xs[2] -= cf->lz;
}

void debug(const char *s) {
	printf("%s",s);
}

void section(const char *s) {
	printf("%s",s);
}

void section_loop(const char *s) {
#ifdef DEBUG_LOOP
	printf("%s",s);
#endif
}

int main(int argc, const char* argv[]) {
	
	struct Config cf = (struct Config){
		.lx = 30,
		.ly = 30,
		.lz = 30,
		.dt = 0.01,
		.tf = 1000000,
		.ndens = 5,
		.mf = 1.0,
		.nbs = 1,
		.rb = 3.0,
		.mb = (1.0*5.0)*(4.0/3.0)*3.141592*3.0*3.0*3.0,
		.Tf = 1.0,
		.fluid_collision_period = 20,
		.output_fluid_config_period = -1,
		.output_big_config_period = 1,
		.output_fluid_velocity_field_measure_period = -1,
		.output_fluid_velocity_field_output_period = -1,
		.fix_big = 0
	};

	int ix, iy, iz, ip, iip, it, ib, ic;

	double lxf = cf.lx;
	double lyf = cf.ly;
	double lzf = cf.lz;
	int nps = cf.ndens*cf.lx*cf.ly*cf.lz;
	sqrttfmf = sqrt(cf.Tf/cf.mf);

	srand(time(NULL));

	// Output things
	FILE * f_config;
	FILE * f_config_big;
	FILE * f_cells;
	FILE * f_vel;
	FILE * f_collisions;
	f_config = fopen("config.txt","w");
	fclose(f_config);
	f_config_big = fopen("config_big.txt","w");
	fclose(f_config_big);
	f_cells = fopen("cells.txt","w");
	fclose(f_cells);
	f_vel = fopen("vel.txt","w");
	fclose(f_vel);
	f_collisions = fopen("collisions.txt","w");
	fclose(f_collisions);

	// Big particles
	section("Allocating big particles\n");
    struct BParticle *bp = (struct BParticle *)malloc(cf.nbs * sizeof(struct BParticle));
    for (ib=0; ib<cf.nbs; ib++) {
		bp[ip].r[0] = cf.lx/2.0;
		bp[ip].r[1] = cf.ly/2.0;
		bp[ip].r[2] = cf.lz/2.0;
		bp[ip].v[0] = 0.0;
		bp[ip].v[1] = 0.0;
		bp[ip].v[2] = 0.0;
		bp[ip].a[0] = 0.0;
		bp[ip].a[1] = 0.0;
		bp[ip].a[2] = 0.0;
		bp[ip].f[0] = 0.0;
		bp[ip].f[1] = 0.0;
		bp[ip].f[2] = 0.0;
		bp[ip].o[0] = 0.0; // particles start pointing in z-direction.
		bp[ip].o[1] = 0.0;
		bp[ip].o[2] = 1.0;
		bp[ip].w[0] = 0.0;
		bp[ip].w[1] = 0.0;
		bp[ip].w[2] = 0.0;
		bp[ip].aa[0] = 0.0;
		bp[ip].aa[1] = 0.0;
		bp[ip].aa[2] = 0.0;
		bp[ip].t[0] = 0.0;
		bp[ip].t[1] = 0.0;
		bp[ip].t[2] = 0.0;
    }
    printf("  .  number of big particles      = %d\n",cf.nbs);

	// Fluid particles
	section("Allocating small particles\n");
	double vt[3] = {0.0,0.0,0.0};
	double vm[3];
	double dist;
	int c = 0;
	int tmp[3];
	int skip = 0;
	struct FParticle *fp = malloc(nps * sizeof(struct FParticle));
	if (fp == NULL) {
		printf("Problema en malloc *fp\n");
		exit(0);
	}
    for (ip=0; ip<nps; ip++) {
    	skip = 0;
		fp[ip].r[0] = randf()*cf.lx;
		fp[ip].r[1] = randf()*cf.ly;
		fp[ip].r[2] = randf()*cf.lz;
		// retry if there is collision with big particle
		for (ib=0; ib<cf.nbs; ib++) {
			dist = distance(fp[ip].r,bp[ib].r,tmp,&cf);
			if (dist < cf.rb) {
				ip -= 1;
				skip = 1;
				break;
			}
		}
		if (skip) continue;
		fp[ip].v[0] = gaussian()*sqrttfmf;
		fp[ip].v[1] = gaussian()*sqrttfmf;
		fp[ip].v[2] = gaussian()*sqrttfmf;
		vt[0] += fp[ip].v[0];
		vt[1] += fp[ip].v[1];
		vt[2] += fp[ip].v[2];
		c += 1;
    }


	printf("  .  number of fluid particles      = %d\n",nps);
	if (nps > 0) {
    	printf("  .  initial total velocity   = %.16f %.16f %.16f\n",vt[0],vt[1],vt[2]);
		printf("  .  initial average velocity = %.16f %.16f %.16f\n",vt[0]/nps,vt[1]/nps,vt[2]/nps);
	}

	section("Substract initial velocity\n");
	// substract initial velocity
	vm[0] = vt[0]/nps;
	vm[1] = vt[1]/nps;
	vm[2] = vt[2]/nps;
	vt[0] = 0.0;
	vt[1] = 0.0;
	vt[2] = 0.0;
	c = 0;
	for (ip=0; ip<nps; ip++) {
		c += 1;
		fp[ip].v[0] -= vm[0];
		fp[ip].v[1] -= vm[1];
		fp[ip].v[2] -= vm[2];
		vt[0] += fp[ip].v[0];
		vt[1] += fp[ip].v[1];
		vt[2] += fp[ip].v[2];
	}
	if (nps > 0) {
		printf("  .  initial total velocity   = %.16f %.16f %.16f\n",vt[0],vt[1],vt[2]);
		printf("  .  initial average velocity = %.16f %.16f %.16f\n",vt[0]/nps,vt[1]/nps,vt[2]/nps);
	}

	// Cells
	section("Allocating cells\n");
	int nmax = 50;
	struct Cell ***cell = (struct Cell ***)malloc(cf.lx * sizeof(struct Cell **));
    for (ix=0; ix<cf.lx; ix++) {
    	cell[ix] = (struct Cell **)malloc(cf.ly * sizeof(struct Cell *));
    	for (iy=0; iy<cf.ly; iy++) {
    		cell[ix][iy] = (struct Cell *)malloc(cf.lz * sizeof(struct Cell));
    		for (iz=0; iz<cf.lz; iz++) {
    			cell[ix][iy][iz].c = 0;
    			cell[ix][iy][iz].t = 0;
    			cell[ix][iy][iz].vc = 0;
    			if (ix < 2 || ix > cf.lx-3 || iy < 2 || iy > cf.ly-3 || iz < 2 || iz > cf.lz-3)
					cell[ix][iy][iz].b = 1;
				else
					cell[ix][iy][iz].b = 0;
    			for (ic=0; ic<nmax; ic++) {
    				cell[ix][iy][iz].p[ic] = -1;
    			}
    			for (ic=0; ic<3; ic++) {
    				cell[ix][iy][iz].v[ic] = 0.0;
    			}
    		}
    	}
    }

    // 
	section("Printing initial configuration\n");
	output_fluid_configuration(fp,f_config,nps);

	section("Starting loop\n");
    // .------.
    // | LOOP |
    // '------'
	double dtt;
	double vmc[DIM] = {0.0,0.0,0.0};
	double vrand[DIM]; // random velocity
	double vrandm[DIM]; // mean random velocity
	double mi[DIM][DIM] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}; // moment of inertia tensor
	double mii[DIM][DIM]; // inverse moment of inertia
	double am[DIM],amr[DIM]; // angular momentum, and cross r
	double rcm[DIM],dv[DIM],dr[DIM]; // center of mass, relative velocity and relative position
	int cc; // cell count
	// rotation things
	double mib, wm[3], wn;
	double pitch, roll, yaw, norm;

	// for timing performance
	double t_cycle = 0.0;
	double t_advection = 0.0;
	clock_t t_advection_begin, t_advection_end;
	double t_advection_collision = 0.0;
	double t_collision = 0.0;
	double t_cells = 0.0;
	clock_t t_cells_begin, t_cells_end;
	double t_output = 0.0
	;

	clock_t t_total_begin = clock();

	for (it=1; it<=cf.tf; it++) {

		clock_t t_cycle_begin = clock();

		if (it%100==0)
			printf("TIMESTEP: %d\n",it);

		section_loop("reset cells\n");
		t_cells_begin = clock();
		for (ix=0; ix<cf.lx; ix++) {
    		for (iy=0; iy<cf.ly; iy++) {
    			for (iz=0; iz<cf.lz; iz++) {
    				cell[ix][iy][iz].c = 0;
    				cell[ix][iy][iz].t = 0;
    				for (ip=0; ip<nmax; ip++) {
    					cell[ix][iy][iz].p[ip] = -1;
    				}
    			}
    		}
    	}
    	t_cells_end = clock();
	    t_cells += (double)(t_cells_end - t_cells_begin)/CLOCKS_PER_SEC;

		section_loop("relevant cells for collision\n");
		t_cells_begin = clock();
    	for (ib=0; ib<cf.nbs; ib++) {
			for (ix=(int)floor(bp[ib].r[0]-cf.rb-2.0); ix<(int)ceil(bp[ib].r[0]+cf.rb+2.0); ix++) {
				for (iy=(int)floor(bp[ib].r[1]-cf.rb-2.0); iy<(int)ceil(bp[ib].r[1]+cf.rb+2.0); iy++) {
					for (iz=(int)floor(bp[ib].r[2]-cf.rb-2.0); iz<(int)ceil(bp[ib].r[2]+cf.rb+2.0); iz++) {
						//printf("%d %d %d %d %d %d\n",ix,iy,iz,mod(ix,lx),mod(iy,ly),mod(iz,lz));
						cell[mod(ix,cf.lx)][mod(iy,cf.ly)][mod(iz,cf.lz)].t = 1;
					}
				}
			}
		}
		t_cells_end = clock();
	    t_cells += (double)(t_cells_end - t_cells_begin)/CLOCKS_PER_SEC;

		// Advection of fluid particles, collision with big ones.
		section_loop("advection fluid particles\n");
		t_advection_begin = clock();
		int shift[3];
        double nc[3], rfc[3], rbc[3], vfn[3];
		for (ip=0; ip<nps; ip++) {

			ix = (int)(fp[ip].r[0]);
			iy = (int)(fp[ip].r[1]);
			iz = (int)(fp[ip].r[2]);

			if (ix < 0 || ix >= cf.lx || iy < 0 || iy >= cf.ly || iz < 0 || iz >= cf.lz) {
				printf("%d %d %d\n",ix,iy,iz);
			}

			fp[ip].r[0] += fp[ip].v[0]*cf.dt;
			fp[ip].r[1] += fp[ip].v[1]*cf.dt;
			fp[ip].r[2] += fp[ip].v[2]*cf.dt;

			if (cell[ix][iy][iz].b == 1) {
				if (fp[ip].r[0] < 0.0)
					fp[ip].r[0] += cf.lx;
				else if (fp[ip].r[0] >= cf.lx)
					fp[ip].r[0] -= cf.lx;
				if (fp[ip].r[1] < 0.0)
					fp[ip].r[1] += cf.ly;
				else if (fp[ip].r[1] >= cf.ly)
					fp[ip].r[1] -= cf.ly;
				if (fp[ip].r[2] < 0.0)
					fp[ip].r[2] += cf.lz;
				else if (fp[ip].r[2] >= cf.lz)
					fp[ip].r[2] -= cf.lz;
			}

			//clock_t t_advection_collision_begin = clock();
        	for (ib=0; ib<cf.nbs; ib++) {
				if (cell[ix][iy][iz].t == 0) continue;
				dist = distance(fp[ip].r,bp[ib].r,shift,&cf);
				if (dist < cf.rb) {
					dtt = collision_position_time(fp[ip].r, bp[ib].r, fp[ip].v, bp[ib].v, shift, rfc, rbc, nc, &cf);
					collision_velocities(rfc, fp[ip].v, rbc, bp[ib].v, bp[ib].w, bp[ib].o, nc, vfn, it, &cf);
					collision_forces(fp[ip].v, vfn, rfc, rbc, bp[ib].f, bp[ib].t, &cf);
					update_srd(dtt, rfc, vfn, fp[ip].r, fp[ip].v, &cf);

				}
			}
			//clock_t t_advection_collision_end = clock();
			//t_advection_collision += (double)(t_advection_collision_end-t_advection_collision_begin)/CLOCKS_PER_SEC;
		}
		t_advection_end = clock();
		t_advection += (double)(t_advection_end-t_advection_begin)/CLOCKS_PER_SEC;

		section_loop("movement of big\n");
		// movement of big particles
		t_advection_begin = clock();
		if (!cf.fix_big) {
			for (ib=0; ib<cf.nbs; ib++) {
				bp[ib].r[0] += bp[ib].v[0]*cf.dt + 0.5*bp[ib].a[0]*cf.dt*cf.dt;
				bp[ib].r[1] += bp[ib].v[1]*cf.dt + 0.5*bp[ib].a[1]*cf.dt*cf.dt;
				bp[ib].r[2] += bp[ib].v[2]*cf.dt + 0.5*bp[ib].a[2]*cf.dt*cf.dt;

				bp[ib].v[0] += 0.5*(bp[ib].a[0] + bp[ib].f[0]/cf.mb)*cf.dt;
				bp[ib].v[1] += 0.5*(bp[ib].a[1] + bp[ib].f[1]/cf.mb)*cf.dt;
				bp[ib].v[2] += 0.5*(bp[ib].a[2] + bp[ib].f[2]/cf.mb)*cf.dt;

				bp[ib].a[0] = bp[ib].f[0]/cf.mb;
				bp[ib].a[1] = bp[ib].f[1]/cf.mb;
				bp[ib].a[2] = bp[ib].f[2]/cf.mb;

				// Angulos
				mib = (2.0/5.0)*cf.mb*cf.rb*cf.rb;

				bp[ib].o[0] += (  bp[ib].o[2]*bp[ib].w[1] - bp[ib].o[1]*bp[ib].w[2])*cf.dt;
				bp[ib].o[1] += (- bp[ib].o[2]*bp[ib].w[0] + bp[ib].o[0]*bp[ib].w[2])*cf.dt;
				bp[ib].o[2] += (  bp[ib].o[1]*bp[ib].w[0] - bp[ib].o[0]*bp[ib].w[1])*cf.dt;

				// ... normalizar vector directors
				norm = sqrt(bp[ib].o[0]*bp[ib].o[0] + bp[ib].o[1]*bp[ib].o[1] + bp[ib].o[2]*bp[ib].o[2]);
				bp[ib].o[0] /= norm;
				bp[ib].o[1] /= norm;
				bp[ib].o[2] /= norm;

    			bp[ib].w[0] += 0.5*(bp[ib].aa[0] + bp[ib].t[0]/mib)*cf.dt;
    			bp[ib].w[1] += 0.5*(bp[ib].aa[1] + bp[ib].t[1]/mib)*cf.dt;
    			bp[ib].w[2] += 0.5*(bp[ib].aa[2] + bp[ib].t[2]/mib)*cf.dt;

				bp[ib].aa[0] = bp[ib].t[0]/mib;
				bp[ib].aa[1] = bp[ib].t[1]/mib;
				bp[ib].aa[2] = bp[ib].t[2]/mib;

				if (bp[ib].r[0] < 0.0)
					bp[ib].r[0] += cf.lx;
				else if (bp[ib].r[0] >= cf.lx)
					bp[ib].r[0] -= cf.lx;
				if (bp[ib].r[1] < 0.0)
					bp[ib].r[1] += cf.ly;
				else if (bp[ib].r[1] >= cf.ly)
					bp[ib].r[1] -= cf.ly;
				if (bp[ib].r[2] < 0.0)
					bp[ib].r[2] += cf.lz;
				else if (bp[ib].r[2] >= cf.lz)
					bp[ib].r[2] -= cf.lz;
			}
		}
		t_advection_end = clock();
		t_advection += (double)(t_advection_end-t_advection_begin)/CLOCKS_PER_SEC;

		section_loop("reset forces\n");
		t_cells_begin = clock();
    	for (ib=0; ib<cf.nbs; ib++) {
    		
			bp[ib].f[0] = 0.0;
			bp[ib].f[1] = 0.0;
			bp[ib].f[2] = 0.0;

			bp[ib].t[0] = 0.0;
			bp[ib].t[1] = 0.0;
			bp[ib].t[2] = 0.0;

			for (ix=(int)floor(bp[ib].r[0]-cf.rb-2.0); ix<(int)ceil(bp[ib].r[0]+cf.rb+2.0); ix++) {
				for (iy=(int)floor(bp[ib].r[1]-cf.rb-2.0); iy<(int)ceil(bp[ib].r[1]+cf.rb+2.0); iy++) {
					for (iz=(int)floor(bp[ib].r[2]-cf.rb-2.0); iz<(int)ceil(bp[ib].r[2]+cf.rb+2.0); iz++) {
						//printf("%d %d %d %d %d %d\n",ix,iy,iz,mod(ix,lx),mod(iy,ly),mod(iz,lz));
						cell[mod(ix,cf.lx)][mod(iy,cf.ly)][mod(iz,cf.lz)].t = 1;
					}
				}
			}
		}
		t_cells_end = clock();
	    t_cells += (double)(t_cells_end - t_cells_begin)/CLOCKS_PER_SEC;


		// -------------------------
		// collision fluid particles
		// -------------------------
		section_loop("collision of fluid particles\n");
		clock_t t_collision_begin = clock();
		if (it%cf.fluid_collision_period == 0) {

			// shift particles
			double cshift[3];
			cshift[0] = randf()-0.5;
			cshift[1] = randf()-0.5;
			cshift[2] = randf()-0.5;
			for (ip=0; ip<nps; ip++) {
				fp[ip].r[0] += cshift[0];
				fp[ip].r[1] += cshift[1];
				fp[ip].r[2] += cshift[2];
				if (fp[ip].r[0] < 0.0)
					fp[ip].r[0] += cf.lx;
				else if (fp[ip].r[0] >= cf.lx)
					fp[ip].r[0] -= cf.lx;
				if (fp[ip].r[1] < 0.0)
					fp[ip].r[1] += cf.ly;
				else if (fp[ip].r[1] >= cf.ly)
					fp[ip].r[1] -= cf.ly;
				if (fp[ip].r[2] < 0.0)
					fp[ip].r[2] += cf.lz;
				else if (fp[ip].r[2] >= cf.lz)
					fp[ip].r[2] -= cf.lz;
				ix = (int)(fp[ip].r[0]);
				iy = (int)(fp[ip].r[1]);
				iz = (int)(fp[ip].r[2]);
				if (cell[ix][iy][iz].c == nmax-1) {
					printf("Out of bounds\n");
					exit(0);
				}
				cell[ix][iy][iz].p[cell[ix][iy][iz].c] = ip;
				cell[ix][iy][iz].c += 1;
			}

			for (ix=0; ix<cf.lx; ix++) {
	    		for (iy=0; iy<cf.ly; iy++) {
	    			for (iz=0; iz<cf.lz; iz++) {

	    				cc = cell[ix][iy][iz].c;

	    				if (cc == 0)
	    					continue;

	    				// mean random velocity
	    				vrandm[0] = 0.0;
	    				vrandm[1] = 0.0;
	    				vrandm[2] = 0.0;

	    				// mean velocity
	    				vmc[0] = 0.0;
	    				vmc[1] = 0.0;
	    				vmc[2] = 0.0;

	    				// center of mass
	    				rcm[0] = 0.0;
	    				rcm[1] = 0.0;
	    				rcm[2] = 0.0;
	    				for (iip=0; iip<cc; iip++) {
	    					ip = cell[ix][iy][iz].p[iip];
	    					rcm[0] += fp[ip].r[0];
	    					rcm[1] += fp[ip].r[1];
	    					rcm[2] += fp[ip].r[2];
	    				}
	    				rcm[0] /= cc;
	    				rcm[1] /= cc;
	    				rcm[2] /= cc;

					    for (iip=0; iip<cc; iip++) {
					    	ip = cell[ix][iy][iz].p[iip];
					    	// ... mean velocity
					    	vmc[0] += fp[ip].v[0];
	    					vmc[1] += fp[ip].v[1];
	    					vmc[2] += fp[ip].v[2];
	    					// ... random velocity
							vrand[0] = gaussian()*sqrttfmf;
							vrand[1] = gaussian()*sqrttfmf;
							vrand[2] = gaussian()*sqrttfmf;
							// ... mean of random velocities
							vrandm[0] += vrand[0];
							vrandm[1] += vrand[1];
							vrandm[2] += vrand[2];
							// ... center of mass
							rcm[0] += fp[ip].r[0];
							rcm[1] += fp[ip].r[1];
							rcm[2] += fp[ip].r[2];
#ifdef AMC
							// ... angular momentum
							//		. relative velocity
							dv[0] = fp[ip].v[0] - vrand[0];
							dv[1] = fp[ip].v[1] - vrand[1];
							dv[2] = fp[ip].v[2] - vrand[2];
							//		. relative position
							dr[0] = fp[ip].r[0] - rcm[0];
							dr[1] = fp[ip].r[1] - rcm[1];
							dr[2] = fp[ip].r[2] - rcm[2];
							//		. cross product between the two
							am[0] = -dr[2]*dv[1] + dr[1]*dv[2];
							am[1] = +dr[2]*dv[0] - dr[0]*dv[2];
							am[2] = -dr[1]*dv[0] + dr[0]*dv[1];
							//		. moment of inertia tensor (with respect to center of mass)
							mi[0][0] += cf.mf*(dr[1]*dr[1] + dr[2]*dr[2]);
							mi[1][1] += cf.mf*(dr[0]*dr[0] + dr[2]*dr[2]);
							mi[2][2] += cf.mf*(dr[0]*dr[0] + dr[1]*dr[1]);
							mi[0][1] -= cf.mf*dr[0]*dr[1];
							mi[0][2] -= cf.mf*dr[0]*dr[2];
							mi[1][0] -= cf.mf*dr[1]*dr[0];
							mi[1][2] -= cf.mf*dr[1]*dr[2];
							mi[2][0] -= cf.mf*dr[2]*dr[0];
							mi[2][1] -= cf.mf*dr[2]*dr[1];
							matrix_inverse(mi,mii);
#endif
							// set velocity to random velocities for now, so we dont have to store them
							fp[ip].v[0] = vrand[0];
							fp[ip].v[1] = vrand[1];
							fp[ip].v[2] = vrand[2];
					    }
					    // means
	    				vmc[0] /= cc;
	    				vmc[1] /= cc;
	    				vmc[2] /= cc;
	    				vrandm[0] /= cc;
	    				vrandm[1] /= cc;
	    				vrandm[2] /= cc;
	    				rcm[0] /= cc;
	    				rcm[1] /= cc;
	    				rcm[2] /= cc;
					    for (iip=0; iip<cc; iip++) {
					    	ip = cell[ix][iy][iz].p[iip];
					    	// final cross product
					    	amr[0] = - am[2]*fp[ip].v[1] + am[1]*fp[ip].v[2];
					    	amr[1] =   am[2]*fp[ip].v[0] - am[0]*fp[ip].v[2];
					    	amr[2] = - am[1]*fp[ip].v[0] + am[0]*fp[ip].v[1];
					    	// aaaand the sheebang
#ifdef AMC
					    	fp[ip].v[0] += vmc[0] - vrandm[0] + cf.mf*mii[0][0]*amr[0] + mii[0][1]*amr[1] + mii[0][2]*amr[2];
					    	fp[ip].v[1] += vmc[1] - vrandm[1] + cf.mf*mii[1][0]*amr[0] + mii[1][1]*amr[1] + mii[1][2]*amr[2];
					    	fp[ip].v[2] += vmc[2] - vrandm[2] + cf.mf*mii[2][0]*amr[0] + mii[2][1]*amr[1] + mii[2][2]*amr[2];
#else
							fp[ip].v[0] += vmc[0] - vrandm[0];
					    	fp[ip].v[1] += vmc[1] - vrandm[1];
					    	fp[ip].v[2] += vmc[2] - vrandm[2];
#endif
					    }
	    			}
	    		}
	    	}

	    	// shift particles back
		    for (ip=0; ip<nps; ip++) {
				fp[ip].r[0] -= cshift[0];
				fp[ip].r[1] -= cshift[1];
				fp[ip].r[2] -= cshift[2];
				if (fp[ip].r[0] < 0.0)
					fp[ip].r[0] += cf.lx;
				else if (fp[ip].r[0] >= cf.lx)
					fp[ip].r[0] -= cf.lx;
				if (fp[ip].r[1] < 0.0)
					fp[ip].r[1] += cf.ly;
				else if (fp[ip].r[1] >= cf.ly)
					fp[ip].r[1] -= cf.ly;
				if (fp[ip].r[2] < 0.0)
					fp[ip].r[2] += cf.lz;
				else if (fp[ip].r[2] >= cf.lz)
					fp[ip].r[2] -= cf.lz;
			}

	    }
	    clock_t t_collision_end = clock();
	    t_collision += (double)(t_collision_end - t_collision_begin)/CLOCKS_PER_SEC;

	    // Output
		section_loop("output\n");
	    clock_t t_output_begin = clock();
	    if (cf.output_fluid_velocity_field_measure_period > 0 && it%cf.output_fluid_velocity_field_measure_period == 0) {
	    	measure_fluid_velocity_field(cell,fp,&cf);
	    }
	    if (cf.output_fluid_velocity_field_output_period > 0 && it%cf.output_fluid_velocity_field_output_period == 0) {
	    	output_fluid_velocity_field(cell,f_vel,&cf);
	    }
		if (cf.output_fluid_config_period > 0 && it%cf.output_fluid_config_period == 0) {
			output_fluid_configuration(fp,f_config,nps);
		}
		if (cf.output_big_config_period > 0 && it%cf.output_big_config_period == 0) {
			output_big_configuration(bp,f_config_big,cf.nbs);
		}
		clock_t t_output_end = clock();
		t_output += (double)(t_output_end - t_output_begin)/CLOCKS_PER_SEC;

		clock_t t_cycle_end = clock();
		t_cycle += (double)(t_cycle_end-t_cycle_begin)/CLOCKS_PER_SEC;
	}

	clock_t t_total_end = clock();
	double t_total = (double)(t_total_end - t_total_begin) / CLOCKS_PER_SEC;
	double t_pcycle = t_cycle/cf.tf;
	printf("Total time: %f\n",t_total);
	printf("Time per cycle: %f\n",t_pcycle);
	printf(" %% advection-collision: %.2f\n",100.0*t_advection/cf.tf/t_pcycle);
	//printf("Time advection collision: %f\n",t_advection_collision/tf);
	printf(" %% rotation: %.2f\n",100.0*t_collision/cf.tf/t_pcycle);
	printf(" %% cells: %.2f\n",100.0*t_cells/cf.tf/t_pcycle);
	printf(" %% output: %.2f\n",100.0*t_output/cf.tf/t_pcycle);
}

	    				// random direction
	    				/*
	    				r1 = randf();
	    				r2 = randf();
	    				phi = 2.0*M_PI*r1;
	    				rho = 2.0*r2-1.0;
	    				rv[0] = cos(phi)*sqrt(1-rho*rho);
	    				rv[1] = sin(phi)*sqrt(1-rho*rho);
	    				rv[2] = rho;
	    				*/

							/*
							// Get relative velocity parallel to random vector
							ur_pn = ur[0]*rv[0]+ur[1]*rv[1]+ur[2]*rv[2];
							ur_pa[0] = ur_pn*rv[0];
							ur_pa[1] = ur_pn*rv[1];
							ur_pa[2] = ur_pn*rv[2];

							// Get relative velocity perpendicular to random vector
							ur_pe[0] = ur[0]-ur_pa[0];
							ur_pe[1] = ur[1]-ur_pa[1];
							ur_pe[2] = ur[2]-ur_pa[2];

							eps[0] = ur_pe[0]*cos(theta) + (rv[2]*ur_pe[1]-rv[1]*ur_pe[2])*sin(theta);
							eps[1] = ur_pe[1]*cos(theta) + (rv[0]*ur_pe[2]-rv[2]*ur_pe[0])*sin(theta);
							eps[2] = ur_pe[2]*cos(theta) + (rv[1]*ur_pe[0]-rv[0]*ur_pe[1])*sin(theta);

							u[0] = eps[0] + ur_pa[0];
							u[1] = eps[1] + ur_pa[1];
							u[2] = eps[2] + ur_pa[2];

							//printf("%f %f\n",rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2],u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
							*/