#include "xy_system.h"

//заполнение массива x, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void fill_first_x(double * first_x, double * first_gvelocity, double * first_dvelocity, ParticleParams params)
{
    for(int i = 0; i < params.amount; ++i)
    {
        first_x[i] = first_gvelocity[i] - first_dvelocity[i];
    }
}

//заполнение массива y, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void fill_first_y(double * first_y, double * first_gvelocity, double * first_dvelocity, double * first_grho,
             double * first_drho, ParticleParams params)
{
    for(int i = 0; i < params.amount; ++i)
    {
        first_y[i] = first_gvelocity[i] + first_drho[i] / first_grho[i] * first_dvelocity[i];
    }
}

void fill_drho_in_xg(double * drho_in_xg, double * x_g, double * image_x_d, double * image_dmass,
                ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    for(int i = 0; i < amount; ++i)
    {
        drho_in_xg[i] = 0;
        for(int j = 0; j < 3 * amount - 2; ++j)
        {
            drho_in_xg[i] += image_dmass[j] * spline_kernel(x_g[i], image_x_d[j], problem_params);
        }
    }
}

double found_next_x(double prev_x, double prev_r, double prev_grho, double prev_drho_xg, double * image_gmass, double * image_dmass,
                    double * prev_image_grho, double * prev_image_r, double * prev_image_x_d,
                    ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double K = problem_params.K;
    double c_s = problem_params.c_s;

    double result = 0;
    double sum = 0;
    double denom = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(prev_r, prev_image_r[j],
                                                                                             problem_params);
    }
    denom = 1. + tau * K * (1. / prev_grho + 1. / prev_drho_xg);
    result = (- tau * pow(c_s, 2) * sum + prev_x) / denom;

    return result;
}

double found_next_y(double prev_y, double prev_r, double prev_grho, double * image_gmass, double * prev_image_grho,
                    double * prev_image_r, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double c_s = problem_params.c_s;

    double result;
    double sum;

    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(prev_r, prev_image_r[j],
                                                                                             problem_params);
    }
    result = - tau * pow(c_s, 2) * sum + prev_y;

    return result;
}

double found_next_dvelocity_xg(double next_x, double next_y, double prev_grho, double prev_drho_xg)
{
    return (next_y - next_x) / (1. + prev_drho_xg / prev_grho);
}

double found_next_gvelocity(double next_x, double next_y, double prev_grho, double prev_drho)
{
    return next_x + found_next_dvelocity_xg(next_x, next_y, prev_grho, prev_drho);
}

void xy_system(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
{
    int gamount = gas_params.amount;
    int damount = dust_params.amount;

    double T = problem_params.T;
    double tau = problem_params.tau;

    //Блок массивов для газа.BEGIN
    double prev_x_g[gamount];
    double prev_gvelocity[gamount];
    double prev_grho[gamount];
    double gmass[gamount];

    double prev_image_x_g[3*gamount - 2];
    double prev_image_gvelocity[3*gamount - 2];
    double prev_image_grho[3*gamount - 2];
    double image_gmass[3*gamount - 2];

    double next_x_g[gamount];
    double next_gvelocity[gamount];
    double next_grho[gamount];
    double next_image_x_g[3 * gamount - 2];
    double next_image_gvelocity[3 * gamount - 2];
    double next_image_ghro[3 * gamount - 2];

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
    double prev_dvelocity[damount];
    double prev_drho[damount];
    double dmass[damount];

    double prev_image_x_d[3 * damount - 2];
    double prev_image_dvelocity[3 * damount - 2];
    double prev_image_drho[3 * damount - 2];
    double image_dmass[3 * damount - 2];

    double next_x_d[damount];
    double next_dvelocity[damount];
    double next_drho[damount];
    double next_image_x_d[3 * damount - 2];
    double next_image_dvelocity[3 * damount - 2];
    double next_image_dhro[3 * damount - 2];

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

    double prev_drho_xg[damount];
    double next_dvelocity_xg[damount];
    fill_drho_in_xg(prev_drho_xg, prev_x_g, prev_image_x_d, image_dmass, gas_params, problem_params);

    //Массивы для x,y.BEGIN
    double prev_x[gamount];
    double prev_y[gamount];

    double next_x[gamount];
    double next_y[gamount];
    
    fill_first_x(prev_x, prev_gvelocity, prev_dvelocity, gas_params);
    fill_first_y(prev_y, prev_gvelocity, prev_dvelocity, prev_grho, prev_drho, gas_params);
    //Массивы для x,y.END

    char fileName[512];
    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
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
            next_x[i] = found_next_x(prev_x[i], prev_x_g[i], prev_grho[i], prev_drho_xg[i], image_gmass, image_dmass,
                                     prev_image_grho, prev_image_x_g, prev_image_x_d, gas_params, problem_params);
            next_y[i] = found_next_y(prev_y[i], prev_x_g[i], prev_grho[i], image_gmass, prev_image_grho, prev_image_x_g,
                                     gas_params, problem_params);

            next_dvelocity_xg[i] = found_next_dvelocity_xg(next_x[i], next_y[i], prev_grho[i], prev_drho_xg[i]);
            next_gvelocity[i] = found_next_gvelocity(next_x[i], next_y[i], prev_grho[i], prev_drho[i]);
            //next_dvelocity[i] = записать сюда аппроксимацию скорости пыли через скорость пыли в точках газа
        }
    }
}