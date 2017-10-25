#include "common_use.h"
#include "gas_stair.h"
#include "gas_wave.h"
#include "dust_wave.h"
#include "four_eq.h"
#include "xy_system.h"
#include "fscanf.h"

int main()
{
    ProblemParams problem_params;
    problem_params.T = 0.5;
    problem_params.h = 0.01;
    problem_params.tau = 0.001;
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

    /*
    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=100_T=0.out",
             "/home/calat/CLionProjects/particles/dustywaveK=1000-t0.000.out");

    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=100_T=0.5.out",
             "/home/calat/CLionProjects/particles/dustywaveK=100-t0.500.out");

    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=1000_T=0.out",
             "/home/calat/CLionProjects/particles/dustywaveK=1000-t0.000.out");

    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=1000_T=0.5.out",
             "/home/calat/CLionProjects/particles/dustywaveK=1000-t0.500.out");

    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=1_T=0.out",
             "/home/calat/CLionProjects/particles/dustywave-t0.000.out");

    xy_value("/home/calat/CLionProjects/particles/xy_value/xy_K=1_T=0.5.out",
             "/home/calat/CLionProjects/particles/dustywave-t0.500.out");

    */

    //stair_gas_print(gs_particle, problem_params);
    //only_dust_wave(dw_particle, problem_params);
    //only_gas_wave(gw_particle, problem_params);
    whole_system(gw_particle, dw_particle, problem_params);

    //xy_system(gw_particle, dw_particle, problem_params);

    return 0;
}
