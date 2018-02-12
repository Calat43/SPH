#pragma once

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "common_use.h"
#include "gas_wave.h"
#include "dust_wave.h"

void near_scheme(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params);

//void smooth(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params);

