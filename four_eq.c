#include "four_eq.h"

double found_next_gvelocity(double x_g, double prev_grho, double prev_gvelocity, double prev_dvelocity,
                            double * image_gmass, double * prev_image_x_g, double * prev_image_grho,
                            ParticleParams gas_params, ProblemParams problem_params)
{
    int amount = gas_params.amount;
    double c_s = problem_params.c_s;
    double tau = problem_params.tau;
    double K = problem_params.K;

    double result = 0;

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_gmass[j] * pow(c_s, 2) * (1. / prev_image_grho[j] + 1. / prev_grho) *
                  spline_gradient(x_g, prev_image_x_g[j], problem_params);
    }
    result = -tau * result;
    result = result - tau * K / prev_grho * (prev_gvelocity - prev_dvelocity) + prev_gvelocity;
    return result;
}

double found_next_dvelocity(double prev_drho, double prev_dvelocity, double prev_gvelocity, ProblemParams params)
{
    double tau = params.tau;
    double K = params.K;

    double result = 0;

    result = tau * K / prev_drho * (prev_gvelocity - prev_dvelocity) + prev_dvelocity;
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

    coordinate_distribution(prev_x_g, gas_params);
    fill_image_x(prev_image_x_g, gas_params);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, gas_params, problem_params);
    fill_image(image_gmass, gmass, gas_params);

    fill_initial_rho(prev_grho, image_gmass, prev_x_g, prev_image_x_g, gas_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, gas_params);

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

    coordinate_distribution(prev_x_d, dust_params);
    fill_image_x(prev_image_x_d, dust_params);

    double average_drho = ddensity_distribution(0);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, dust_params, problem_params);
    fill_image(image_dmass, dmass, dust_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params);

    fill_image(prev_image_drho, prev_drho, dust_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, dust_params);
    //Блок массивов для пыли.END

    double dvelocity;
    double gvelocity;
    char fileName[512];
    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        dvelocity = 0;
        gvelocity = 0;

        sprintf(fileName, "/home/calat/CLionProjects/particles/grho/grho_%d.dat", frameId);
        FILE * grho_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(grho_frame, "%lf %0.15lf\n", prev_x_g[i], prev_grho[i]);
        }
        fclose(grho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/gvelocity/gvel_%d.dat", frameId);
        FILE * gvelocity_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(gvelocity_frame, "%lf %0.15lf\n", prev_x_g[i], prev_gvelocity[i]);
        }
        fclose(gvelocity_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/drho/drho_%d.dat", frameId);
        FILE * drho_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(drho_frame, "%lf %0.15lf\n", prev_x_d[i], prev_drho[i]);
        }
        fclose(drho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/dvelocity/dvel_%d.dat", frameId);
        FILE * dvelocity_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(dvelocity_frame, "%lf %0.15lf\n", prev_x_d[i], prev_dvelocity[i]);
        }
        fclose(dvelocity_frame);

        for(int i = 0; i < gamount; ++i)
        {
            dvelocity = point_value(prev_x_g[i], prev_image_dvelocity, image_dmass, prev_image_drho, prev_image_x_d,
                                    dust_params, problem_params);
            next_grho[i] = found_next_rho(image_gmass, prev_x_g, prev_image_x_g, i, gas_params, problem_params);
            next_gvelocity[i] = found_next_gvelocity(prev_x_g[i], prev_grho[i], prev_gvelocity[i], dvelocity, image_gmass,
                                                     prev_image_x_g, prev_image_grho, gas_params, problem_params);
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }

        for(int i = 0; i < damount; ++i)
        {
            gvelocity = point_value(prev_x_d[i], prev_image_gvelocity, image_gmass, prev_image_grho, prev_image_x_g,
                                    gas_params, problem_params);
            next_drho[i] = found_next_rho(image_dmass, prev_x_d, prev_image_x_d, i, dust_params, problem_params);
            next_dvelocity[i] = found_next_dvelocity(prev_drho[i], prev_dvelocity[i], gvelocity, problem_params);
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }

        fill_image(next_image_grho, next_grho, gas_params);
        fill_image(next_image_gvelocity, next_gvelocity, gas_params);

        fill_image(next_image_drho, next_drho, dust_params);
        fill_image(next_image_dvelocity, next_dvelocity, dust_params);

        for (int i = 0; i < 3 * gamount - 2; ++i)
        {
            next_image_x_g[i] = found_next_coordinate(prev_image_x_g[i], prev_image_gvelocity[i], problem_params);
        }

        for (int i = 0; i < 3 * damount - 2; ++i)
        {
            next_image_x_d[i] = found_next_coordinate(prev_image_x_d[i], prev_image_dvelocity[i], problem_params);
        }

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

}