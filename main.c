#include "common_use.h"
#include "gas_stair.h"
#include "gas_wave.h"
#include "dust_wave.h"
#include "four_eq.h"
#include "xy_system.h"

int main()
{
    ProblemParams problem_params;
    problem_params.T = 0.5;
    problem_params.h = 0.04;
    problem_params.tau = 0.02;
    problem_params.c_s = 1;
    problem_params.K = 1;

    ParticleParams gs_particle;
    gs_particle.amount = 400;
    gs_particle.left = 0;
    gs_particle.right = 1;

    ParticleParams dw_particle;
    dw_particle.amount = 400;
    dw_particle.left = 0;
    dw_particle.right = 1;

    ParticleParams gw_particle;
    gw_particle.amount = 400;
    gw_particle.left = 0;
    gw_particle.right = 1;

    //stair_gas_print(gs_particle, problem_params);
    //only_dust_wave(dw_particle, problem_params);
    //only_gas_wave(gw_particle, problem_params);
    //whole_system(gw_particle, dw_particle, problem_params);

    xy_system(gw_particle, dw_particle, problem_params);

    return 0;
}
