#pragma once

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "common_use.h"
#include "gas_wave.h"
#include "dust_wave.h"

void near(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params);

void smooth(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params);

