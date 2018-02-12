#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include "common_use.h"
#include "gas_stair.h"
#include "gas_wave.h"
#include "dust_wave.h"
#include "explicit.h"
#include "xy_system.h"
#include "fscanf.h"
#include "cells.h"

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

    //explicit_scheme(gw_particle, dw_particle, problem_params);
    //near_scheme(gw_particle, dw_particle, problem_params);
    cells(problem_params.h, gw_particle, dw_particle, problem_params);
    

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
    problem_params.delta = 1. / 10000;
    problem_params.t_stop = problem_params.d2g / problem_params.K;


    //ATTENTION! (dust amount = gas amount) || (dust amount = 2 * gas amount) || (2 * dust amount = gas amount)
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

    // Create DATA_DIR if not exists
    //DATA_DIR and PROBLEM_PARAMS_FILE change in common_use.c
    struct stat st;
    if (stat(DATA_DIR, &st) == -1) {
        mkdir(DATA_DIR, 0700);
    }

    FILE * paramsFile = fopen(PROBLEM_PARAMS_FILE, "r");
    if (paramsFile == NULL) {
        printf(stderr, "Error opening file: %s\n", strerror( errno ));
        return errno;
    }
    while (true)
    {
        int params_read = fscanf(
                paramsFile,
                "%lf %lf %lf %lf %lf",
                &(problem_params.delta),
                &(problem_params.d2g),
                &(problem_params.K),
                &(problem_params.h),
                &(problem_params.tau)
        );
        if (params_read != 5) {
            break;
        }

        problem_params.t_stop = problem_params.d2g / problem_params.K;

        solve_problem(dw_particle, gw_particle, problem_params);
   }


    return 0;
}
