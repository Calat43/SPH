#include <time.h>
#include <assert.h>
#include "common_use.h"
#include "gas_stair.h"
#include "gas_wave.h"
#include "dust_wave.h"
#include "four_eq.h"
#include "xy_system.h"
#include "fscanf.h"

void solve_problem(ParticleParams dw_particle, ParticleParams gw_particle, ProblemParams problem_params)
{
    printf(
            "Using problem parameters: d2g: %lf K: %lf h: %lf tau: %lf\n",
            problem_params.d2g,
            problem_params.K,
            problem_params.h,
            problem_params.tau
    );

    clock_t startTime = clock();

    ParticleParams gs_particle;
    gs_particle.amount = 400;
    gs_particle.left = 0;
    gs_particle.right = 1;

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
    //whole_system(gw_particle, dw_particle, problem_params);

    near(gw_particle, dw_particle, problem_params);
    //smooth(gw_particle, dw_particle, problem_params);

    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);
}

int main()
{
    ProblemParams problem_params;
    problem_params.T = 0.5;
    problem_params.h = 0.05;
    problem_params.tau = 0.005;
    problem_params.c_s = 1;
    problem_params.K = 5;
    problem_params.d2g = 0.01;
    problem_params.middle_rho_gas = 1;
    problem_params.delta = 1. / 100;

    //параметры пыли
    ParticleParams dw_particle;
    dw_particle.amount = 400;
    dw_particle.left = 0;
    dw_particle.right = 1;
    dw_particle.isGas = false;

    //параметры газа
    ParticleParams gw_particle;
    gw_particle.amount = 200;
    gw_particle.left = 0;
    gw_particle.right = 1;
    gw_particle.isGas = true;


    FILE * paramsFile = fopen("problem_params.txt", "r");
    while (true)
    {
        int params_read = fscanf(
                paramsFile,
                "%lf %lf %lf %lf",
                &(problem_params.d2g),
                &(problem_params.K),
                &(problem_params.h),
                &(problem_params.tau)
        );
        if (params_read != 4) {
            break;
        }

        solve_problem(dw_particle, gw_particle, problem_params);
   }


    return 0;
}
