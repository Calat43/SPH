#pragma once

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "common_use.h"

void only_dust_wave(ParticleParams particle_params, ProblemParams problem_params);

//начальное распределение плотности пыли
double ddensity_distribution(double x);

//заполнение массива, содержащего массы частиц пыли
void fill_dmass(double * dmass, double * x_d, double * image_x_d, double average_drho,
                ParticleParams particle_params, ProblemParams problem_params);

void fill_initial_dvelocity(double * dvelocity, double * x_d, ParticleParams params);