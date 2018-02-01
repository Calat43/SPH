#pragma once

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

extern const double pi;

extern const char * DATA_DIR;
extern const char * PROBLEM_PARAMS_FILE;


typedef struct _particleParams{
    int amount;
    double left;
    double right;
    bool isGas;
} ParticleParams;

typedef struct _problemParams{
    double T;
    double h;
    double tau;
    double c_s;
    double K;
    double t_stop;
    double d2g;
    double middle_rho_gas;
    double delta;
} ProblemParams;

double spline_kernel(double x_a, double x_b, ProblemParams params);

double spline_gradient(double x_a, double x_b, ProblemParams params);

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params);

//void coordinate_distribution(double * x_d, ParticleParams params);

//void fill_image_x(double * image_x, ParticleParams params);

void fill_x(double * coord, ProblemParams problemParams, ParticleParams particleParams);

void fill_image_x(double * image_coord, double * coord, ParticleParams params);

void fill_image(double * image, double * real, ParticleParams params);

void fill_mass(double * mass, ParticleParams params, ProblemParams problemParams);

void fill_initial_rho(double * rho, double  * image_mass, double * x, double * image_x,
                      ParticleParams particle_params, ProblemParams problem_params);

double found_next_rho(double * image_mass, double * x, double * image_x, int i,
                            ParticleParams particle_params, ProblemParams problem_params);

//x - точка, в которой ищем значение
double interpolation_value(double x, double * image_function, double * image_mass, double * image_rho, double * image_x,
                    ParticleParams particle_params, ProblemParams problem_params);

double image_interpolation_value(double * image_x, double * image_function, double * image_mass, double * image_rho,
                                 double * interpol_x, int i, ParticleParams particleParams, ProblemParams problemParams);

double interpolation_value_for_rho(double x, double * image_mass, double * image_x, int i, ParticleParams particle_params,
                           ProblemParams problem_params);