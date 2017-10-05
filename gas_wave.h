#pragma once

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "common_use.h"

void only_gas_wave(ParticleParams particle_params, ProblemParams problem_params);

//начальное распределение плотности газа
double gdensity_distribution(double x);

//заполнение массива, содержащего массы частиц газа
void fill_gmass(double * gmass, double * x_g, double * image_x_g, double average_grho,
                ParticleParams particle_params, ProblemParams problem_params);

void fill_initial_gvelocity(double * gvelocity, double * x_g, ParticleParams params);
