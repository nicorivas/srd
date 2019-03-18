#pragma once

struct Config {
	int lx;
	int ly;
	int lz;
	double dt;
	int tf;
	int ndens;
	double mf;
	int nbs;
	double rb;	// the big particle radius
	double mb;
	double Tf;
	int fluid_collision_period;
	int output_fluid_config_period;
	int output_big_config_period;
	int output_fluid_velocity_field_measure_period;
	int output_fluid_velocity_field_output_period;
	int fix_big;
};