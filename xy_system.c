#include "xy_system.h"
#include <assert.h>
#include <stdbool.h>

double dust_vel(double x, double t, double d2g)
{
    return d2g * sin((x + t)/2.*pi);
}

//заполнение массива x, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void compute_x(double * x, double * gvelocity, double * dvelocity, ParticleParams params)
{
    for(int i = 0; i < params.amount; ++i)
    {
        x[i] = gvelocity[i] - dvelocity[i];
    }
}

//заполнение массива y, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void compute_y(double * y, double * gvelocity, double * dvelocity, double * grho, double * drho, ParticleParams params,
               ProblemParams problemParams)
{
    for(int i = 0; i < params.amount; ++i)
    {
        y[i] = gvelocity[i] + problemParams.d2g * dvelocity[i];
    }
}

//x^{n+1}, r - координата искомой точки
double found_next_x(double prev_x, double * image_gmass, double * prev_image_grho, double * prev_image_x_g,
                    double prev_grho, double prev_drho, double r, int i, ParticleParams particle_params,
                    ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double K = problem_params.K;
    double c_s = problem_params.c_s;
    double d2g = problem_params.d2g;
    double t_stop = problem_params.t_stop;
    int hN = (int)floor(problem_params.h * amount);

    double result = 0;
    double sum = 0;
    double denom = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(r, prev_image_x_g[j],
                                                                                             problem_params);
    }
    denom = 1. + tau * (d2g / t_stop + 1. / t_stop);
    result = (- tau * pow(c_s, 2) * sum + prev_x) / denom;

    return result;
}

//y^{n+1}, r - координата искомой точки
double found_next_y(double prev_y, double * image_gmass, double * prev_image_grho, double * prev_image_x_g,
                    double prev_grho, double r, int i, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double c_s = problem_params.c_s;
    int hN = (int)floor(problem_params.h * amount);

    double result = 0;
    double sum = 0;

    for(int j = 0; j < 3 * amount - 2; ++j)
    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(r, prev_image_x_g[j],
                                                                                             problem_params);
    }
    result = - tau * pow(c_s, 2) * sum + prev_y;

    return result;
}

static double found_next_dvelocity(double next_x, double next_y, ProblemParams params)
{
    double d2g = params.d2g;
    return (next_y - next_x) / (1. + d2g);
}

static double found_next_gvelocity(double next_x, double next_y, ProblemParams params)
{
    double d2g = params.d2g;
    return (next_y + d2g * next_x)  / (1 + d2g);
}

double near_velocity(double r, double * coordinate, double * velocity, ParticleParams params)
{
    int nearly = 0;


    for (int i = 0; i < params.amount; ++i)
    {
        double new_fabs = fabs(r - coordinate[i]);
        double old_fabs = fabs(r - coordinate[nearly]);
        if (new_fabs < old_fabs)
        {
            nearly = i;
        }
    }

    return velocity[nearly];
}

//осторожно с границами массивов
//значение для газа в точках пыли
double near_for_gas(double r, double * coordinate, int i_gas, double * dust_function,
                    double gamount, double damount, ParticleParams params)
{
    if(gamount == damount)
    {
        return near_velocity(r, coordinate, dust_function, params);
    }
    if(2 *gamount == damount)
    {
        return (dust_function[2 * i_gas] + dust_function[2 * i_gas + 1]) / 2;
    }
    if(gamount == 2 * damount)
    {
        if(i_gas % 2 == 0)
        {
            return dust_function[i_gas / 2];
        }
        if(i_gas % 2 != 0)
        {
            return dust_function[(i_gas - 1) / 2];
        }
        assert(false);
    }
    assert(false);
}

//значение для пыли в точках газа
double near_for_dust(double r, double * coordinate, int i_dust, double * gas_function,
                     double gamount, double damount, ParticleParams params)
{
    if(gamount == damount)
    {
        return near_velocity(r, coordinate, gas_function, params);
    }
    if(2 * gamount == damount)
    {
        if (i_dust % 2 == 0)
        {
            return gas_function[i_dust / 2];
        }
        if (i_dust % 2 != 0)
        {
            return gas_function[(i_dust - 1) / 2];
        }
        assert(false);
    }
    if(gamount == 2 * damount)
    {
        return(gas_function[2 * i_dust] + gas_function[2 * i_dust + 1]) / 2;
    }
    assert(false);
}

void near(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
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

    double prev_image_x_g[3 * gamount - 2];
    double prev_image_gvelocity[3 * gamount - 2];
    double prev_image_grho[3 * gamount - 2];
    double image_gmass[3 * gamount - 2];

    double next_x_g[gamount];
    double next_gvelocity[gamount];
    double next_grho[gamount];
    double next_image_x_g[3 * gamount - 2];
    double next_image_gvelocity[3 * gamount - 2];
    double next_image_grho[3 * gamount - 2];

    for (int i = 0; i < gamount; ++i)
    {
        prev_x_g[i] = NAN;
        prev_gvelocity[i] = NAN;
        prev_grho[i] = NAN;
        gmass[i] = NAN;
        next_x_g[i] = NAN;
        next_gvelocity[i] = NAN;
        next_grho[i] = NAN;
    }

    for (int j = 0; j < 3 * gamount - 2; ++j)
    {
        prev_image_x_g[j] = NAN;
        prev_image_gvelocity[j] = NAN;
        prev_image_grho[j] = NAN;
        image_gmass[j] = NAN;
        next_image_x_g[j] = NAN;
        next_image_gvelocity[j] = NAN;
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
    double next_image_drho[3 * damount - 2];

    for (int i = 0; i < damount; ++i)
    {
        prev_x_d[i] = NAN;
        prev_dvelocity[i] = NAN;
        prev_drho[i] = NAN;
        dmass[i] = NAN;
        next_x_d[i] = NAN;
        next_dvelocity[i] = NAN;
        next_drho[i] = NAN;
    }

    for(int j = 0; j < 3 * damount - 2; ++j)
    {
        prev_image_x_d[j] = NAN;
        prev_image_dvelocity[j] = NAN;
        prev_image_drho[j] = NAN;
        image_dmass[j] = NAN;
        next_image_x_d[j] = NAN;
        next_image_dvelocity[j] = NAN;
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

    //плотность  пыли, посчитанная в точках газа
    double prev_drho_xg[gamount];
    double next_drho_xg[gamount];
    //плотность газа, посчитанная в точка пыли
    double prev_grho_xd[damount];
    double next_grho_xd[damount];

    double prev_gvel_xd[damount];
    double prev_dvel_xg[gamount];

    for(int i = 0; i < gamount; ++i)
    {
        prev_drho_xg[i] = NAN;
        next_drho_xg[i] = NAN;
        prev_dvel_xg[i] = NAN;
    }
    for(int i = 0; i < damount; ++i)
    {
        prev_grho_xd[i] = NAN;
        next_grho_xd[i] = NAN;
        prev_gvel_xd[i] = NAN;
    }

    for(int i = 0; i < damount; ++i)
    {
        //prev_grho_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_grho, gas_params);
        prev_grho_xd[i] = near_for_dust(prev_x_d[i], prev_x_g, i, prev_grho, gamount, damount, gas_params);
        //prev_gvel_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_gvelocity, gas_params);
        prev_gvel_xd[i] = near_for_dust(prev_x_d[i], prev_x_g, i, prev_gvelocity, gamount, damount, gas_params);
    }
    for(int i = 0; i < gamount; ++i)
    {
        //prev_drho_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_drho, dust_params);
        prev_drho_xg[i] = near_for_gas(prev_x_g[i], prev_x_d, i, prev_drho, gamount, damount, dust_params);
        //prev_dvel_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_dvelocity, dust_params);
        prev_dvel_xg[i] = near_for_gas(prev_x_g[i], prev_x_d, i, prev_dvelocity, gamount, damount, dust_params);
    }

    //Массивы для x,y.BEGIN
    double prev_x_gas[gamount];
    double prev_y_gas[gamount];

    double prev_x_dust[damount];
    double prev_y_dust[damount];

    double next_x_gas[gamount];
    double next_y_gas[gamount];

    double next_x_dust[damount];
    double next_y_dust[damount];

    for(int i = 0; i < damount; ++i)
    {
        prev_x_dust[i] = NAN;
        prev_y_dust[i] = NAN;
        next_x_dust[i] = NAN;
        next_y_dust[i] = NAN;
    }
    for(int i = 0; i < gamount; ++i)
    {
        prev_x_gas[i] = NAN;
        prev_y_gas[i] = NAN;
        next_x_gas[i] = NAN;
        next_y_gas[i] = NAN;
    }

    compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
    compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params, problem_params);

    compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
    compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params, problem_params);
    //Массивы для x,y.END

    for (int i = 0; i < gamount; ++i)
    {
        assert(!isnan(prev_x_g[i]));
        assert(!isnan(prev_gvelocity[i]));
        assert(!isnan(prev_grho[i]));
        assert(!isnan(gmass[i]));
    }

    for (int j = 0; j < 3 * gamount - 2; ++j)
    {
        assert(!isnan(prev_image_x_g[j]));
        assert(!isnan(prev_image_gvelocity[j]));
        assert(!isnan(prev_image_grho[j]));
        assert(!isnan(image_gmass[j]));
    }

    for (int i = 0; i < damount; ++i)
    {
        assert(!isnan(prev_x_d[i]));
        assert(!isnan(prev_dvelocity[i]));
        assert(!isnan(prev_drho[i]));
        assert(!isnan(dmass[i]));
    }

    for(int j = 0; j < 3 * damount - 2; ++j)
    {
        assert(!isnan(prev_image_x_d[j]));
        assert(!isnan(prev_image_dvelocity[j]));
        assert(!isnan(prev_image_drho[j]));
        assert(!isnan(image_dmass[j]));
    }

    for(int i = 0; i < gamount; ++i)
    {
        assert(!isnan(prev_drho_xg[i]));
        assert(!isnan(prev_dvel_xg[i]));
    }
    for(int i = 0; i < damount; ++i)
    {
        assert(!isnan(prev_grho_xd[i]));
        assert(!isnan(prev_gvel_xd[i]));
    }

    for(int i = 0; i < damount; ++i)
    {
        assert(!isnan(prev_x_dust[i]));
        assert(!isnan(prev_y_dust[i]));
    }
    for(int i = 0; i < gamount; ++i)
    {
        assert(!isnan(prev_x_gas[i]));
        assert(!isnan(prev_y_gas[i]));
    }

    char fileName[512];

    sprintf(fileName, "%s/xy_gas_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_gas_0 = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_0, "%lf %lf %lf\n", prev_x_g[i], prev_gvelocity[i], prev_grho[i]);
    }
    fclose(xy_gas_0);

    sprintf(fileName, "%s/xy_dust_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_dust_0 = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(xy_dust_0, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(xy_dust_0);

    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        printf("%d\n", frameId);

/*
        sprintf(fileName, "%s/gas/xy_gas_%d_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
                frameId, problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
        FILE * xy_gas_frame = fopen(fileName, "w");
        for (int j = 0; j < 3 * gamount - 2; ++j)
        {
            fprintf(xy_gas_frame, "%lf %lf %lf\n", prev_image_x_g[j], prev_image_gvelocity[j], prev_image_grho[j]);
        }
        fclose(xy_gas_frame);

        sprintf(fileName, "%s/dust/xy_dust_%d_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
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
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
            assert(!isnan(next_x_g[i]));
        }
        for(int i = 0; i < damount; ++i)
        {
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
            assert(!isnan(next_x_d[i]));
        }

        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            next_image_x_g[j] = found_next_coordinate(prev_image_x_g[j], prev_image_gvelocity[j], problem_params);
            assert(!isnan(next_image_x_g[j]));
        }
        for(int j = 0; j < 3 * damount - 2; ++j)
        {
            next_image_x_d[j] = found_next_coordinate(prev_image_x_d[j], prev_image_dvelocity[j], problem_params);
            assert(!isnan(next_image_x_d[j]));
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_grho[i] = found_next_rho(image_gmass, next_x_g, next_image_x_g, i, gas_params, problem_params);
            assert(!isnan(next_grho[i]));
        }
        for(int i = 0; i < damount; ++i)
        {
            next_drho[i] = found_next_rho(image_dmass, next_x_d, next_image_x_d, i, dust_params, problem_params);
            assert(!isnan(next_drho[i]));
        }

        for(int i = 0; i < gamount; ++i)
        {
            //next_drho_xg[i] = near_velocity(next_x_g[i], next_x_d, next_drho, dust_params);
            next_drho_xg[i] = near_for_gas(next_x_g[i], next_x_d, i, next_drho, gamount, damount, dust_params);
            assert(!isnan(next_drho_xg[i]));
        }
        for(int i = 0; i < damount; ++i)
        {
            //next_grho_xd[i] = near_velocity(next_x_d[i], next_x_g, next_grho, gas_params);
            next_grho_xd[i] = near_for_dust(next_x_d[i], next_x_g, i, next_grho, gamount, damount, gas_params);
            assert(!isnan(next_grho_xd[i]));
        }

        //ищем x, y в точках газа
        for(int i = 0; i < gamount; ++i)
        {
            next_x_gas[i] = found_next_x(prev_x_gas[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho[i],
                                     prev_drho_xg[i], prev_x_g[i], i, gas_params, problem_params);
            assert(!isnan(next_x_gas[i]));
            next_y_gas[i] = found_next_y(prev_y_gas[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho[i],
                                     prev_x_g[i], i, gas_params, problem_params);
            assert(!isnan(next_y_gas[i]));
            next_gvelocity[i] = found_next_gvelocity(next_x_gas[i], next_y_gas[i], problem_params);
            assert(!isnan(next_gvelocity[i]));
        }

        //ищем x, y  в точках пыли
        for(int i = 0; i < damount; ++i)
        {
            next_x_dust[i] = found_next_x(prev_x_dust[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho_xd[i],
                                     prev_drho[i], prev_x_d[i], i, gas_params, problem_params);
            assert(!isnan(next_x_dust[i]));
            next_y_dust[i] = found_next_y(prev_y_dust[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho_xd[i],
                                     prev_x_d[i], i, gas_params, problem_params);
            assert(!isnan(next_y_dust[i]));
            next_dvelocity[i] = found_next_dvelocity(next_x_dust[i], next_y_dust[i], problem_params);
            assert(!isnan(next_dvelocity[i]));

        }

        fill_image(next_image_drho, next_drho, dust_params);
        fill_image(next_image_grho, next_grho, gas_params);
        fill_image(next_image_dvelocity, next_dvelocity, dust_params);
        fill_image(next_image_gvelocity, next_gvelocity, gas_params);

        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            assert(!isnan(next_image_gvelocity[j]));
            assert(!isnan(next_image_grho[j]));
        }

        for(int j = 0; j < 3 * damount - 2; ++j)
        {
            assert(!isnan(next_image_dvelocity[j]));
            assert(!isnan(next_image_drho[j]));
        }

        for(int i = 0; i < gamount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];

            prev_x_gas[i] = next_x_gas[i];
            prev_y_gas[i] = next_y_gas[i];

            prev_x_dust[i] = next_x_dust[i];
            prev_y_dust[i] = next_y_dust[i];

            prev_drho_xg[i] = next_drho_xg[i];
        }
        for(int i = 0; i < damount; ++i)
        {
            prev_drho[i] = next_drho[i];
            prev_dvelocity[i] = next_dvelocity[i];
            prev_x_d[i] = next_x_d[i];

            prev_x_dust[i] = next_x_dust[i];
            prev_y_dust[i] = next_y_dust[i];

            prev_grho_xd[i] = next_grho_xd[i];
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

        for(int i = 0; i < gamount; ++i)
        {
            //prev_dvel_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_dvelocity, gas_params);
            prev_dvel_xg[i] = near_for_gas(prev_x_g[i], prev_x_d, i, prev_dvelocity, gamount, damount, gas_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            //prev_gvel_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_gvelocity, dust_params);
            prev_gvel_xd[i] = near_for_dust(prev_x_d[i], prev_x_g, i, prev_gvelocity, gamount, damount, dust_params);
        }

        compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
        compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params, problem_params);

        compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
        compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params, problem_params);
    }

    sprintf(fileName, "%s/xy_gas_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_gas_T = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_T, "%lf %lf %lf\n",prev_x_g[i], prev_gvelocity[i],prev_grho[i]);
    }
    fclose(xy_gas_T);

    sprintf(fileName, "%s/xy_dust_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_dust_T = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(xy_dust_T, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(xy_dust_T);
}

/*
void smooth(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
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

    double prev_image_x_g[3 * gamount - 2];
    double prev_image_gvelocity[3 * gamount - 2];
    double prev_image_grho[3 * gamount - 2];
    double image_gmass[3 * gamount - 2];

    double next_x_g[gamount];
    double next_gvelocity[gamount];
    double next_grho[gamount];
    double next_image_x_g[3 * gamount - 2];
    double next_image_gvelocity[3 * gamount - 2];
    double next_image_grho[3 * gamount - 2];

    for (int i = 0; i < gamount; ++i)
    {
        prev_x_g[i] = NAN;
        prev_gvelocity[i] = NAN;
        prev_grho[i] = NAN;
        gmass[i] = NAN;
        next_x_g[i] = NAN;
        next_gvelocity[i] = NAN;
        next_grho[i] = NAN;
    }

    for (int j = 0; j < 3 * gamount - 2; ++j)
    {
        prev_image_x_g[j] = NAN;
        prev_image_gvelocity[j] = NAN;
        prev_image_grho[j] = NAN;
        image_gmass[j] = NAN;
        next_image_x_g[j] = NAN;
        next_image_gvelocity[j] = NAN;
        next_image_grho[j] = NAN;
    }

    //coordinate_distribution(prev_x_g, gas_params);
    //fill_image_x(prev_image_x_g, gas_params);
    fill_x(prev_x_g, problem_params, gas_params);
    fill_image_x(prev_image_x_g, prev_x_g, gas_params);

    double average_grho = gdensity_distribution(0, problem_params);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, 0, gas_params, problem_params);
    fill_image(image_gmass, gmass, gas_params);

    fill_initial_rho(prev_grho, image_gmass, prev_x_g, prev_image_x_g, gas_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, gas_params, problem_params);

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
    double next_image_drho[3 * damount - 2];

    for (int i = 0; i < damount; ++i)
    {
        prev_x_d[i] = NAN;
        prev_dvelocity[i] = NAN;
        prev_drho[i] = NAN;
        dmass[i] = NAN;
        next_x_d[i] = NAN;
        next_dvelocity[i] = NAN;
        next_drho[i] = NAN;
    }

    for(int j = 0; j < 3 * damount - 2; ++j)
    {
        prev_image_x_d[j] = NAN;
        prev_image_dvelocity[j] = NAN;
        prev_image_drho[j] = NAN;
        image_dmass[j] = NAN;
        next_image_x_d[j] = NAN;
        next_image_dvelocity[j] = NAN;
        next_image_drho[j] = NAN;
    }

    //coordinate_distribution(prev_x_d, dust_params);
    //fill_image_x(prev_image_x_d, dust_params);
    fill_x(prev_x_d, problem_params, dust_params);
    fill_image_x(prev_image_x_d, prev_x_d, dust_params);

    double average_drho = ddensity_distribution(0, problem_params);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, 0, dust_params, problem_params);
    fill_image(image_dmass, dmass, dust_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params, problem_params);

    fill_image(prev_image_drho, prev_drho, dust_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, dust_params);
    //Блок массивов для пыли.END

    //плотность  пыли, посчитанная в точках газа
    double prev_drho_xg[damount];
    double next_drho_xg[damount];
    //плотность газа, посчитанная в точка пыли
    double prev_grho_xd[gamount];
    double next_grho_xd[gamount];

    double prev_gvel_xd[gamount];
    double prev_dvel_xg[gamount];

    for(int i = 0; i < damount; ++i)
    {
        prev_grho_xd[i] = NAN;
        next_grho_xd[i] = NAN;
        prev_gvel_xd[i] = NAN;
    }
    for(int i = 0; i < gamount; ++i)
    {
        prev_drho_xg[i] = NAN;
        next_drho_xg[i] = NAN;
        prev_dvel_xg[i] = NAN;
    }

    for(int i = 0; i < gamount; ++i)
    {
        prev_drho_xg[i] = interpolation_value_for_rho(prev_x_g[i], image_dmass, prev_image_x_d, i, dust_params,
                                                      problem_params);

        //prev_gvel_xd[i] = dust_vel(prev_x_d[i], 0, problem_params.d2g);
        prev_dvel_xg[i] = interpolation_value(prev_x_g, prev_image_dvelocity, image_dmass, prev_image_drho,
                                              prev_image_x_d, i, dust_params, problem_params);
    }
    for(int i = 0; i < damount; ++i)
    {
        prev_grho_xd[i] = interpolation_value_for_rho(prev_x_d[i], image_gmass, prev_image_x_g, i, gas_params,
                                                      problem_params);
        prev_gvel_xd[i] = interpolation_value(prev_x_d, prev_image_gvelocity, image_gmass, prev_image_grho,
                                              prev_image_x_g, i, gas_params, problem_params);
    }

    //Массивы для x,y.BEGIN
    double prev_x_gas[gamount];
    double prev_y_gas[gamount];

    double prev_x_dust[damount];
    double prev_y_dust[damount];

    double next_x_gas[gamount];
    double next_y_gas[gamount];

    double next_x_dust[damount];
    double next_y_dust[damount];

    for(int i = 0; i < damount; ++i)
    {
        prev_x_dust[i] = NAN;
        prev_y_dust[i] = NAN;
        next_x_dust[i] = NAN;
        next_y_dust[i] = NAN;
    }

    for(int i = 0; i < gamount; ++i)
    {
        prev_x_gas[i] = NAN;
        prev_y_gas[i] = NAN;
        next_x_gas[i] = NAN;
        next_y_gas[i] = NAN;
    }

    compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
    compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params);

    compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
    compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params);
    //Массивы для x,y.END

    char fileName[512];

    sprintf(fileName, "%s/sm_xy_gas_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
    FILE * xy_gas_0 = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_0, "%lf %lf %lf\n", prev_x_g[i], prev_gvelocity[i], prev_grho[i]);
    }
    fclose(xy_gas_0);

    sprintf(fileName, "%s/sm_xy_dust_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * xy_dust_0 = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(xy_dust_0, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(xy_dust_0);

    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        printf("%d\n", frameId);

        for(int i = 0; i < gamount; ++i)
        {
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }

        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            next_image_x_g[j] = found_next_coordinate(prev_image_x_g[j], prev_image_gvelocity[j], problem_params);
        }
        for(int j = 0; j < 3 * damount - 2; ++j)
        {
            next_image_x_d[j] = found_next_coordinate(prev_image_x_d[j], prev_image_dvelocity[j], problem_params);
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_grho[i] = found_next_rho(image_gmass, next_x_g, next_image_x_g, i, gas_params, problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            next_drho[i] = found_next_rho(image_dmass, next_x_d, next_image_x_d, i, dust_params, problem_params);
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_drho_xg[i] = interpolation_value_for_rho(next_x_g[i], image_dmass, next_image_x_d, i,
                                                          dust_params, problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            next_grho_xd[i] = interpolation_value_for_rho(next_x_d[i], image_gmass, next_image_x_g, i,
                                                          gas_params, problem_params);
        }

        //ищем x, y в точках газа
        for(int i = 0; i < gamount; ++i)
        {
            next_x_gas[i] = found_next_x(prev_x_gas[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho[i],
                                         prev_drho_xg[i], prev_x_g[i], i, gas_params, problem_params);
            next_y_gas[i] = found_next_y(prev_y_gas[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho[i],
                                         prev_x_g[i], i, gas_params, problem_params);
            next_gvelocity[i] = found_next_gvelocity(next_x_gas[i], next_y_gas[i], next_drho_xg[i]/next_grho[i]);
        }

        //ищем x, y  в точках пыли
        for(int i = 0; i < damount; ++i)
        {
            next_x_dust[i] = found_next_x(prev_x_dust[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho_xd[i],
                                          prev_drho[i], prev_x_d[i], i, gas_params, problem_params);
            next_y_dust[i] = found_next_y(prev_y_dust[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho_xd[i],
                                          prev_x_d[i], i, gas_params, problem_params);
            next_dvelocity[i] = found_next_dvelocity(next_x_dust[i], next_y_dust[i], next_drho[i]/next_grho_xd[i]);

        }

        fill_image(next_image_drho, next_drho, dust_params);
        fill_image(next_image_grho, next_grho, gas_params);
        fill_image(next_image_dvelocity, next_dvelocity, dust_params);
        fill_image(next_image_gvelocity, next_gvelocity, gas_params);

        for(int i = 0; i < gamount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];

            prev_x_gas[i] = next_x_gas[i];
            prev_y_gas[i] = next_y_gas[i];

            prev_drho_xg[i] = next_drho_xg[i];
        }
        for(int i = 0; i < damount; ++i)
        {
            prev_drho[i] = next_drho[i];
            prev_dvelocity[i] = next_dvelocity[i];
            prev_x_d[i] = next_x_d[i];

            prev_x_dust[i] = next_x_dust[i];
            prev_y_dust[i] = next_y_dust[i];

            prev_grho_xd[i] = next_grho_xd[i];
        }

        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            prev_image_grho[j] = next_image_grho[j];
            prev_image_gvelocity[j] = next_image_gvelocity[j];
            prev_image_x_g[j] = next_image_x_g[j];
        }
        for(int j = 0; j < 3 * damount - 2; ++j)
        {
            prev_image_drho[j] = next_image_drho[j];
            prev_image_dvelocity[j] = next_image_dvelocity[j];
            prev_image_x_d[j] = next_image_x_d[j];
        }

        for(int i = 0; i < gamount; ++i)
        {
            //prev_gvel_xd[i] = dust_vel(prev_x_d[i], tau * frameId, problem_params.d2g);

            prev_dvel_xg[i] = interpolation_value(prev_x_g, prev_image_dvelocity, image_dmass, prev_image_drho,
                                                  prev_image_x_d, i, dust_params, problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            prev_gvel_xd[i] = interpolation_value(prev_x_d, prev_image_gvelocity, image_gmass, prev_image_grho,
                                                  prev_image_x_g, i, gas_params, problem_params);
        }

        compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
        compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params);

        compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
        compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params);
    }

    sprintf(fileName, "%s/sm_xy_gas_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ngas=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, gas_params.amount);
    FILE * xy_gas_T = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_T, "%lf %lf %lf\n",prev_x_g[i], prev_gvelocity[i],prev_grho[i]);
    }
    fclose(xy_gas_T);

    sprintf(fileName, "%s/sm_xy_dust_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf_Ndust=%d.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K, dust_params.amount);
    FILE * xy_dust_T = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(xy_dust_T, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(xy_dust_T);
}
*/