#pragma once

#include <stdio.h>
#include <math.h>

const double pi = 3.14159265358;

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