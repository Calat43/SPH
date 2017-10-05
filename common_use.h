#pragma once

#include <stdio.h>
#include <math.h>

extern const double pi;

typedef struct _particleParams{
    int amount;
    double left;
    double right;
} ParticleParams;

typedef struct _problemParams{
    double T;
    double h;
    double tau;
    double c_s;
} ProblemParams;

double spline_kernel(double x_a, double x_b, ProblemParams params);

double spline_gradient(double x_a, double x_b, ProblemParams params);

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params);

void coordinate_distribution(double * x_d, ParticleParams params);

double fill_image_x(double * image_x, ParticleParams params);

void fill_image(double * image, double * real, ParticleParams params);

void fill_initial_rho(double * rho, double  * image_mass, double * x, double * image_x,
                      ParticleParams particle_params, ProblemParams problem_params);

double found_next_image_rho(double * image_mass, double * x, double * image_x, int i,
                            ParticleParams particle_params, ProblemParams problem_params);