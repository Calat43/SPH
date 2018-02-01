#include <assert.h>
#include "four_eq.h"

static double found_next_gvelocity(double x_g, double prev_grho, double prev_gvelocity, double prev_dvelocity,
                            double * image_gmass, double * prev_image_x_g, double * prev_image_grho, int i,
                            ParticleParams gas_params, ProblemParams problem_params)
{
    int amount = gas_params.amount;
    double c_s = problem_params.c_s;
    double tau = problem_params.tau;
    double K = problem_params.K;
    double d2g = problem_params.d2g;
    double t_stop = problem_params.t_stop;
    int hN = (int)floor(problem_params.h * amount);

    double result = 0;

    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) *
                  spline_gradient(x_g, prev_image_x_g[j], problem_params);
    }
    result = -tau * result * c_s * c_s;
    result = result - tau * d2g / t_stop * (prev_gvelocity - prev_dvelocity) + prev_gvelocity;
    return result;
}

static double found_next_dvelocity(double prev_drho, double prev_dvelocity, double prev_gvelocity, ProblemParams params)
{
    double tau = params.tau;
    double K = params.K;
    double t_stop = params.t_stop;

    double result = 0;

    //result = tau * K / prev_drho * (prev_gvelocity - prev_dvelocity) + prev_dvelocity;
    result = tau / t_stop * (prev_gvelocity - prev_dvelocity) + prev_dvelocity;
    return result;
}

double next_rho(double prev_x, double * prev_image_x, double prev_rho, double * prev_image_vel, double * image_mass,
                ParticleParams particleParams, ProblemParams problemParams)
{
    int amount = particleParams.amount;
    double tau = problemParams.tau;
    double result = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_mass[j] * prev_image_vel[j] * spline_kernel(prev_x, prev_image_x[j], problemParams);
    }
    result = - result * tau * prev_rho + prev_rho;
    return result;
}

void whole_system(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
{
    int gamount = gas_params.amount;
    int damount = dust_params.amount;

    double T = problem_params.T;
    double tau = problem_params.tau;

    //Блок массивов для газа.BEGIN
    double prev_x_g[gamount];
    double next_x_g[gamount];

    double gmass[gamount];

    double prev_gvelocity[gamount];
    double next_gvelocity[gamount];
    double prev_grho[gamount];
    double next_grho[gamount];

    double prev_image_x_g[3*gamount - 2];
    double next_image_x_g[3*gamount - 2];

    double image_gmass[3*gamount - 2];

    double prev_image_gvelocity[3*gamount - 2];
    double next_image_gvelocity[3*gamount - 2];
    double prev_image_grho[3*gamount - 2];
    double next_image_grho[3*gamount - 2];

    for(int i = 0; i < gamount; ++i)
    {
        prev_x_g[i] = NAN;
        next_x_g[i] = NAN;
        gmass[i] = NAN;
        prev_gvelocity[i] = NAN;
        next_gvelocity[i] = NAN;
        prev_grho[i] = NAN;
        next_grho[i] = NAN;
    }

    for(int j = 0; j < 3 * gamount - 2; ++j)
    {
        prev_image_x_g[j] = NAN;
        next_image_x_g[j] = NAN;
        image_gmass[j] = NAN;
        prev_image_gvelocity[j] = NAN;
        next_image_gvelocity[j] = NAN;
        prev_image_grho[j] = NAN;
        next_image_grho[j] = NAN;
    }

    //coordinate_distribution(prev_x_g, gas_params);
    //fill_image_x(prev_image_x_g, gas_params);

    fill_x(prev_x_g, problem_params, gas_params);
    fill_image_x(prev_image_x_g, prev_x_g, gas_params);

    double average_grho = gdensity_distribution(0, problem_params);
    //fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, 0, gas_params, problem_params);
    fill_mass(gmass, gas_params, problem_params);
    fill_image(image_gmass, gmass, gas_params);

    fill_initial_rho(prev_grho, image_gmass, prev_x_g, prev_image_x_g, gas_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, gas_params, problem_params);

    fill_image(prev_image_grho, prev_grho, gas_params);
    fill_image(prev_image_gvelocity, prev_gvelocity, gas_params);
    //Блок массивов для газа.END

    //Блок массивов для пыли.BEGIN
    double prev_x_d[damount];
    double next_x_d[damount];

    double dmass[damount];

    double prev_dvelocity[damount];
    double next_dvelocity[damount];
    double prev_drho[damount];
    double next_drho[damount];

    double prev_image_x_d[3 * damount - 2];
    double next_image_x_d[3 * damount - 2];

    double image_dmass[3 * damount - 2];

    double prev_image_dvelocity[3 * damount - 2];
    double next_image_dvelocity[3 * damount - 2];
    double prev_image_drho[3 * damount - 2];
    double next_image_drho[3 * damount - 2];


    for(int i = 0; i < damount; ++i)
    {
        prev_x_d[i] = NAN;
        next_x_d[i] = NAN;
        dmass[i] = NAN;
        prev_dvelocity[i] = NAN;
        next_dvelocity[i] = NAN;
        prev_drho[i] = NAN;
        next_drho[i] = NAN;
    }

    for(int j = 0; j < 3 * damount - 2; ++j)
    {
        prev_image_x_d[j] = NAN;
        next_image_x_d[j] = NAN;
        image_dmass[j] = NAN;
        prev_image_dvelocity[j] = NAN;
        next_image_dvelocity[j] = NAN;
        prev_image_drho[j] = NAN;
        next_image_drho[j] = NAN;
    }

    //coordinate_distribution(prev_x_d, dust_params);
    //fill_image_x(prev_image_x_d, dust_params);

    fill_x(prev_x_d, problem_params, dust_params);
    fill_image_x(prev_image_x_d, prev_x_d, dust_params);

    double average_drho = ddensity_distribution(0, problem_params);
    //fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, 0, dust_params, problem_params);
    fill_mass(dmass, dust_params, problem_params);
    fill_image(image_dmass, dmass, dust_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params, problem_params);

    fill_image(prev_image_drho, prev_drho, dust_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, dust_params);
    //Блок массивов для пыли.END

    double prev_dvel_xg[gamount];
    double prev_gvel_xd[damount];

    for(int i  = 0; i < gamount; ++i)
    {
        prev_dvel_xg[i] = NAN;
    }
    for(int i = 0; i < damount; ++i)
    {
        prev_gvel_xd[i] = NAN;
    }

    char fileName[512];

    /*
    sprintf(fileName, "%s/explicit_gas_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
    FILE * explicit_gas_0 = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(explicit_gas_0, "%lf %lf %lf\n", prev_x_g[i], prev_gvelocity[i], prev_grho[i]);
    }
    fclose(explicit_gas_0);

    sprintf(fileName, "%s/explicit_dust_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * explicit_dust_0 = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(explicit_dust_0, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(explicit_dust_0);
*/

    sprintf(fileName, "%s/explicit_gas_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);

    FILE * explicit_imgas_0 = fopen(fileName, "w");
    for (int j = 0; j < 3 * gamount - 2; ++j)
    {
        fprintf(explicit_imgas_0, "%lf %lf %lf\n", prev_image_x_g[j], prev_image_gvelocity[j], prev_image_grho[j]);
    }
    fclose(explicit_imgas_0);

    sprintf(fileName, "%s/explicit_dust_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * explicit_imdust_0 = fopen(fileName, "w");
    for (int j = 0; j < 3 * damount - 2; ++j)
    {
        fprintf(explicit_imdust_0, "%lf %lf %lf\n", prev_image_x_d[j], prev_image_dvelocity[j], prev_image_drho[j]);
    }
    fclose(explicit_imdust_0);

    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        printf("%d\n", frameId);

        /*
        sprintf(fileName, "%s/gas/expl_gas_%d_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
                frameId, problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
        FILE * xy_gas_frame = fopen(fileName, "w");
        for (int j = 0; j < 3 * gamount - 2; ++j)
        {
            fprintf(xy_gas_frame, "%lf %lf %lf\n", prev_image_x_g[j], prev_image_gvelocity[j], prev_image_grho[j]);
        }
        fclose(xy_gas_frame);

        sprintf(fileName, "%s/dust/expl_dust_%d_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
                frameId, problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
        FILE * xy_dust_frame = fopen(fileName, "w");
        for (int j = 0; j < 3 * damount - 2; ++j)
        {
            fprintf(xy_dust_frame, "%lf %lf %lf\n", prev_image_x_d[j], prev_image_dvelocity[j], prev_image_drho[j]);
        }
        fclose(xy_dust_frame);
*/

        for(int i = 0; i < gamount; ++i)
        {
            prev_dvel_xg[i] = interpolation_value(prev_x_g[i], prev_image_dvelocity, image_dmass, prev_image_drho,
                                                  prev_image_x_d, dust_params, problem_params);
            assert(!isnan(prev_dvel_xg[i]));
        }
        for(int i = 0; i < damount; ++i)
        {
            prev_gvel_xd[i] = interpolation_value(prev_x_d[i], prev_image_gvelocity, image_gmass, prev_image_grho,
                                                  prev_image_x_g, gas_params, problem_params);
            assert(!isnan(prev_gvel_xd[i]));
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_gvelocity[i] = found_next_gvelocity(prev_x_g[i], prev_grho[i], prev_gvelocity[i], prev_dvel_xg[i], image_gmass,
                                                     prev_image_x_g, prev_image_grho, i, gas_params, problem_params);
            assert(!isnan(next_gvelocity[i]));
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }

        for(int i = 0; i < damount; ++i)
        {
            next_dvelocity[i] = found_next_dvelocity(prev_drho[i], prev_dvelocity[i], prev_gvel_xd[i], problem_params);
            assert(!isnan(next_dvelocity[i]));
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }

        fill_image(next_image_gvelocity, next_gvelocity, gas_params);

        fill_image(next_image_dvelocity, next_dvelocity, dust_params);

        for (int i = 0; i < 3 * gamount - 2; ++i)
        {
            next_image_x_g[i] = found_next_coordinate(prev_image_x_g[i], prev_image_gvelocity[i], problem_params);
            assert(!isnan(next_image_x_g[i]));
        }

        for (int i = 0; i < 3 * damount - 2; ++i)
        {
            next_image_x_d[i] = found_next_coordinate(prev_image_x_d[i], prev_image_dvelocity[i], problem_params);
            assert(!isnan(next_image_x_d[i]));
        }
        for(int i = 0; i < gamount; ++i)
        {
            //next_grho[i] = next_rho(prev_x_g[i], prev_image_x_g, prev_grho[i], prev_image_gvelocity, image_gmass, gas_params, problem_params);
            next_grho[i] = found_next_rho(image_gmass, next_x_g, next_image_x_g, i, gas_params, problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            //next_drho[i] = next_rho(prev_x_d[i], prev_image_x_d, prev_drho[i], prev_image_dvelocity, image_dmass, dust_params, problem_params);
            next_drho[i] = found_next_rho(image_dmass, next_x_d, next_image_x_d, i, dust_params, problem_params);
        }
        fill_image(next_image_grho, next_grho, gas_params);
        fill_image(next_image_drho, next_drho, dust_params);

        for(int i = 0; i < gamount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];
        }

        for(int i = 0; i < damount; ++i)
        {
            prev_drho[i] = next_drho[i];
            prev_dvelocity[i] = next_dvelocity[i];
            prev_x_d[i] = next_x_d[i];
        }

        for(int i = 0; i < 3 * gamount - 2; ++i)
        {
            prev_image_grho[i] = next_image_grho[i];
            prev_image_gvelocity[i] = next_image_gvelocity[i];
            prev_image_x_g[i] = next_image_x_g[i];

        }

        for(int i = 0; i < 3 * damount - 2; ++i)
        {
            prev_image_drho[i] = next_image_drho[i];
            prev_image_dvelocity[i] = next_image_dvelocity[i];
            prev_image_x_d[i] = next_image_x_d[i];

        }
    }

    /*
    sprintf(fileName, "%s/explicit_gas_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
    FILE * explicit_gas_T = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(explicit_gas_T, "%lf %lf %lf\n", prev_x_g[i], prev_gvelocity[i], prev_grho[i]);
    }
    fclose(explicit_gas_T);
    sprintf(fileName, "%s/explicit_dust_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * explicit_dust_T = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(explicit_dust_T, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(explicit_dust_T);
*/
    sprintf(fileName, "%s/explicit_gas_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
    FILE * explicit_imgas_T = fopen(fileName, "w");
    for (int j = 0; j < 3 * gamount - 2; ++j)
    {
        fprintf(explicit_imgas_T, "%lf %lf %lf\n", prev_image_x_g[j], prev_image_gvelocity[j], prev_image_grho[j]);
    }
    fclose(explicit_imgas_T);
    sprintf(fileName, "%s/explicit_dust_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * explicit_imdust_T = fopen(fileName, "w");
    for (int j = 0; j < 3 * damount - 2; ++j)
    {
        fprintf(explicit_imdust_T, "%lf %lf %lf\n", prev_image_x_d[j], prev_image_dvelocity[j], prev_image_drho[j]);
    }
    fclose(explicit_imdust_T);

    printf("----");

}