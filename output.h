#pragma once
#include <stdio.h>
#include "particles.h"
#include "cells.h"
#include "config.h"

void output_fluid_configuration(struct FParticle *fp, FILE *f, int n);
void measure_fluid_velocity_field(struct Cell ***cell, struct FParticle *fp, struct Config *cf);
void output_fluid_velocity_field(struct Cell ***cell, FILE *f, struct Config *cf);
void output_big_configuration(struct BParticle *bp, FILE *f, int n);
void output_cells(struct Cell ***c, FILE *f, struct Config *cf);