#include "output.h"
#include "particles.h"
#include "cells.h"

void output_fluid_configuration(struct FParticle *fp, FILE *f, int n) {
	int ip;
	f = fopen("config.txt","a");
	for (ip=0; ip<n; ip++) {
		if (ip == 10)
			fprintf(f,"%.16f %.16f %.16f %.16f %.16f %.16f\n",fp[ip].r[0],fp[ip].r[1],fp[ip].r[2],fp[ip].v[0],fp[ip].v[1],fp[ip].v[2]);
	}
	fclose(f);
}

void measure_fluid_velocity_field(struct Cell ***cell, struct FParticle *fp, struct Config *cf) {
	int ix,iy,iz,cc,iip,ip;
	double *v;
	for (ix=0; ix<cf->lx; ix++) {
		for (iy=0; iy<cf->ly; iy++) {
			for (iz=0; iz<cf->lz; iz++) {
				cc = cell[ix][iy][iz].c;
				if (cc == 0) continue;
				//v = cell[ix][iy][iz].v;
			    for (iip=0; iip<cc; iip++) {
			    	ip = cell[ix][iy][iz].p[iip];
			    	cell[ix][iy][iz].v[0] += fp[ip].v[0];
					cell[ix][iy][iz].v[1] += fp[ip].v[1];
					cell[ix][iy][iz].v[2] += fp[ip].v[2];
				    cell[ix][iy][iz].vc += 1;
			    }
			}
		}
	}
}

void output_fluid_velocity_field(struct Cell ***cell, FILE *f, struct Config *cf) {
	int ix,iy,iz;
	double v[3];
	f = fopen("vel.txt","a");
	for (ix=0; ix<cf->lx; ix++) {
		for (iy=0; iy<cf->ly; iy++) {
			for (iz=0; iz<cf->lz; iz++) {
				if (cell[ix][iy][iz].vc == 0) {
					fprintf(f,"%d %d %d %.16f %.16f %.16f\n",ix,iy,iz,0.0,0.0,0.0);
					continue;
				}
				v[0] = cell[ix][iy][iz].v[0]/cell[ix][iy][iz].vc;
				v[1] = cell[ix][iy][iz].v[1]/cell[ix][iy][iz].vc;
				v[2] = cell[ix][iy][iz].v[2]/cell[ix][iy][iz].vc;
				fprintf(f,"%d %d %d %.16f %.16f %.16f\n",ix,iy,iz,v[0],v[1],v[2]);
			}
		}
	}
	fclose(f);
}

void output_big_configuration(struct BParticle *bp, FILE *f, int n) {
	int ib;
	f = fopen("config_big.txt","a");
    for (ib=0; ib<n; ib++) {
    	fprintf(f,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n",
    		bp[ib].r[0],bp[ib].r[1],bp[ib].r[2],
    		bp[ib].v[0],bp[ib].v[1],bp[ib].v[2],
    		bp[ib].f[0],bp[ib].f[1],bp[ib].f[2],
    		bp[ib].o[0],bp[ib].o[1],bp[ib].o[2],
    		bp[ib].w[0],bp[ib].w[1],bp[ib].w[2],
    		bp[ib].t[0],bp[ib].t[1],bp[ib].t[2]);
    }
	fclose(f);
}

void output_cells(struct Cell ***c, FILE *f, struct Config *cf) {
	int ix, iy, iz;
	f = fopen("cells.txt","a");
	for (ix=0; ix<cf->lx; ix++) {
		for (iy=0; iy<cf->ly; iy++) {
			for (iz=0; iz<cf->lz; iz++) {
				fprintf(f,"%d %d %d %d\n",ix,iy,iz,c[ix][iy][iz].t);
			}
		}
	}
	fclose(f);
}