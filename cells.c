#include "cells.h"

void fill_cells_num(int * cells_num, double * x, double cell_length, int cell_amount, ParticleParams params)
{
    for (int i = 0; i < params.amount; ++i)
    {
        for(int j = 1; j <= cell_amount; ++j)
        {
            if((x[i] >= params.left + cell_length * (j - 1)) && (x[i] < params.left + cell_length * j))
            {
                cells_num[i] = j - 1;
            }
        }
    }
}

//в первом столбце - количество частиц в ячейке
void fill_cells(int ** particle_cells, int * cell_num, int cell_amount, ParticleParams params)
{
    int j = 0;

    for(int l = 0; l < cell_amount; ++l)
    {
        j = 1;
        for(int i = 0; i < params.amount; ++i)
        {
            if(cell_num[i] == l)
            {
                particle_cells[l][j] = i;
                ++j;
            }
        }
        particle_cells[l][j] = -1;
        particle_cells[l][0] = j;
    }
}

//int cell_num - номер v*
double vel_asterisk(int cell_num, double * velocity, int ** cells)
{
    double result = 0;
    double neighbors = cells[cell_num][0];
    int neighbor = 0;

    for (int k = 1; k <= neighbors; ++k)
    {
        neighbor = cells[cell_num][k];
        result += velocity[neighbor];
    }
    return result / neighbors;
}

double eps_asterisk(int cell_num, double gmass, double dmass, int ** gas_cells, int ** dust_cells)
{
    double gas_neighbors = (double) gas_cells[cell_num][0];
    double dust_neighbors = (double) dust_cells[cell_num][0];

    return (dmass * dust_neighbors) / (gmass * gas_neighbors);
}

double t_stop_asterisk(int cell_num, double K, double * drho, int ** dust_cells)
{
    return K / vel_asterisk(cell_num, drho, dust_cells);
}

double a_p(double grho, double x_g, double gmass, int gamount, double * image_grho, double * image_x_g,
           ProblemParams params)
{
    double result = 0;
    double c_s = params.c_s;

    for(int j = 0; j < 3 * gamount - 2; ++j)
    {
        result += (1. / image_grho[j] + 1. / grho) * spline_gradient(x_g, image_x_g[j], params);
    }

    return - result * gmass * pow(c_s, 2);
}

double a_p_asterisk(int cell_num, int ** cells, double * grho, double * x_g, double gmass, int gamount, double * image_grho,
                double * image_x_g, ProblemParams params)
{
    double neighbors = cells[cell_num][0];
    int neighbor = 0;
    double result = 0;

    for(int k = 1; k <= neighbors; ++k)
    {
        neighbor = cells[cell_num][k];
        result += a_p(grho[neighbor], x_g[neighbor], gmass, gamount, image_grho, image_x_g, params);
    }
    return result / neighbors;
}

static void compute_x(double * x, double * gvel_astr, double * dvel_astr, int cell_amount)
{
    for(int i = 0; i < cell_amount; ++i)
    {
        x[i] = gvel_astr[i] - dvel_astr[i];
    }
}

static void compute_y(double * y, double * gvel_astr, double * dvel_astr, double eps_astr, int cell_amount)
{
    for(int i = 0; i < cell_amount; ++i)
    {
        y[i] = gvel_astr[i] + eps_astr * dvel_astr[i];
    }
}

static double found_next_x(double prev_x, double a_p_astr, double t_stop_astr, double eps_astr, ProblemParams params)
{
    double denom = params.tau / t_stop_astr * (eps_astr - 1);

    return (params.tau * a_p_astr + prev_x) / denom;
}

static double found_next_y(double prev_y, double a_p_astr, ProblemParams params)
{
    return params.tau * a_p_astr + prev_y;
}

double found_next_gvel_astr(double x, double y, double eps_astr)
{
    return (y + eps_astr * x) / (1 + eps_astr);
}

double found_next_dvel_astr(double x, double y, double eps_astr)
{
    return (y - x) / (1 + eps_astr);
}

double next_gvel(int i, double eps_astr, double t_stop_astr, double prev_gvel, double next_dvel_astr, double * grho,
                 double * x_g, double gmass, int gamount, double * image_grho, double * image_x_g, ProblemParams params)
{
    double result = 0;
    double denom = 0;
    double tau = params.tau;
    double ap = a_p(grho[i], x_g[i], gmass, gamount, image_grho, image_x_g, params);

    denom = 1. / tau + eps_astr / t_stop_astr;
    result = (prev_gvel / tau + ap + eps_astr / t_stop_astr * next_dvel_astr) / denom;

    return result;
}

double next_dvel(double prev_dvel, double t_stop_astr, double next_gvel_astr, ProblemParams params)
{
    double result = 0;
    double denom = 0;
    double tau = params.tau;

    denom =  1./ tau + 1 / t_stop_astr;
    result = (prev_dvel / tau + next_gvel_astr / t_stop_astr) / denom;

    return result;
}

void cells(double cell_length, ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
{
    int gamount = gas_params.amount;
    int damount = dust_params.amount;

    double T = problem_params.T;
    double tau = problem_params.tau;

    double gmass = compute_mass(gas_params, problem_params);
    double dmass = compute_mass(dust_params, problem_params);
    int cell_num = 0;

    //Если выражение возвращает не целый тип?
    int cell_amount = (int)((gas_params.right - gas_params.left) / cell_length);

    double a_p_astr = 0;
    double eps_astr = 0;
    double t_stop_astr = 0;

    //arrays for gas. BEGIN
    double prev_x_g[gamount];
    double prev_gvelocity[gamount];
    double prev_grho[gamount];

    double prev_image_x_g[3 * gamount - 2];
    double prev_image_gvelocity[3 * gamount - 2];
    double prev_image_grho[3 * gamount - 2];

    double next_x_g[gamount];
    double next_gvelocity[gamount];
    double next_grho[gamount];

    double next_image_x_g[3 * gamount - 2];
    double next_image_gvelocity[3 * gamount - 2];
    double next_image_grho[3 * gamount - 2];

    fill_x(prev_x_g, problem_params, gas_params);
    fill_image_x(prev_image_x_g, prev_x_g, gas_params);

    fill_initial_rho_const_mass(prev_grho, gmass, prev_x_g, prev_image_x_g, gas_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, gas_params, problem_params);

    fill_image(prev_image_grho, prev_grho, gas_params);
    fill_image(prev_image_gvelocity, prev_gvelocity, gas_params);
    //arrays for gas. END

    //arrays for dust. BEGIN
    double prev_x_d[damount];
    double prev_dvelocity[damount];
    double prev_drho[damount];

    double prev_image_x_d[3 * damount - 2];
    double prev_image_dvelocity[3 * damount - 2];
    double prev_image_drho[3 * damount - 2];

    double next_x_d[damount];
    double next_dvelocity[damount];
    double next_drho[damount];

    double next_image_x_d[3 * damount - 2];
    double next_image_dvelocity[3 * damount - 2];
    double next_image_drho[3 * damount - 2];

    fill_x(prev_x_d, problem_params, dust_params);
    fill_image_x(prev_image_x_d, prev_x_d, dust_params);

    fill_initial_rho_const_mass(prev_drho, dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params, problem_params);

    fill_image(prev_image_drho, prev_drho, dust_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, dust_params);
    //arrays for dust. END

    //Если gamount, damount - нечетное?
    int gas_cells_num[gamount];
    int _gas_cells[cell_amount][gamount / 2];
    int * gas_cells[cell_amount];
    for (int i = 0; i < cell_amount; ++i)
    {
        gas_cells[i] = _gas_cells[i];
    }

    int dust_cells_num[damount];
    int _dust_cells[cell_amount][damount / 2];
    int * dust_cells[cell_amount];
    for (int i = 0; i < cell_amount; ++i)
    {
        dust_cells[i] = _dust_cells[i];
    }

    double prev_gvel_astr[cell_amount];
    double next_gvel_astr[cell_amount];

    double prev_dvel_astr[cell_amount];
    double next_dvel_astr[cell_amount];

    double prev_x[cell_amount];
    double next_x[cell_amount];

    double prev_y[cell_amount];
    double next_y[cell_amount];

    char fileName[512];

    sprintf(fileName, "%s/cell_gas_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_gas_0 = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_0, "%lf %lf %lf\n", prev_x_g[i], prev_gvelocity[i], prev_grho[i]);
    }
    fclose(xy_gas_0);

    sprintf(fileName, "%s/cell_dust_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
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

        fill_cells_num(gas_cells_num, prev_x_g, cell_length, cell_amount, gas_params);
        fill_cells(gas_cells, gas_cells_num, cell_amount, gas_params);

        fill_cells_num(dust_cells_num, prev_x_d, cell_length, cell_amount, dust_params);
        fill_cells(dust_cells, dust_cells_num, cell_amount, dust_params);

        a_p_astr = a_p_asterisk(cell_num, gas_cells, prev_grho, prev_x_g, gmass, gamount, prev_image_grho,
                                       prev_image_x_g, problem_params);
        eps_astr = eps_asterisk(cell_num, gmass, dmass, gas_cells, dust_cells);
        t_stop_astr = t_stop_asterisk(cell_num, problem_params.K, prev_drho, dust_cells);


        for(int k = 0; k < cell_amount; ++k)
        {
            prev_gvel_astr[k] = vel_asterisk(k, prev_gvelocity, gas_cells);
            prev_dvel_astr[k] = vel_asterisk(k, prev_dvelocity, dust_cells);
        }

        compute_x(prev_x, prev_gvel_astr, prev_dvel_astr, cell_amount);
        compute_y(prev_y, prev_gvel_astr, prev_dvel_astr, eps_astr, cell_amount);

        for(int k = 0; k < cell_amount; ++k)
        {
            next_x[k] = found_next_x(prev_x[k], a_p_astr, t_stop_astr, eps_astr, problem_params);
            next_y[k] = found_next_y(prev_y[k], a_p_astr, problem_params);
        }

        for(int k = 0; k < cell_amount; ++k)
        {
            next_dvel_astr[k] = found_next_dvel_astr(next_x[k], next_y[k], eps_astr);
            next_gvel_astr[k] = found_next_gvel_astr(next_x[k], next_y[k], eps_astr);
        }

        for(int i = 0; i < gamount; ++i)
        {
            cell_num = gas_cells_num[i];
            next_gvelocity[i] = next_gvel(i, eps_astr, t_stop_astr, prev_gvelocity[i], next_dvel_astr[cell_num],
                                          prev_grho, prev_x_g, gmass, gamount, prev_image_grho, prev_image_x_g,
                                          problem_params);
       }

        for(int i = 0; i < damount; ++i)
        {
            cell_num = dust_cells_num[i];
            next_dvelocity[i] = next_dvel(prev_dvelocity[i], t_stop_astr, next_gvel_astr[cell_num], problem_params);
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_x_g[i] = next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            next_x_d[i] = next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }
        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            next_image_x_g[j] = next_coordinate(prev_image_x_g[j], prev_image_gvelocity[j], problem_params);
        }
        for(int j = 0; j < 3 * damount - 2; ++j)
        {
            next_image_x_d[j] = next_coordinate(prev_image_x_d[j], prev_image_dvelocity[j], problem_params);
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_grho[i] = next_rho_const_mass(gmass, next_x_g, next_image_x_g, i, gas_params, problem_params);
        }
        for(int i = 0; i < damount; ++i)
        {
            next_drho[i] = next_rho_const_mass(dmass, next_x_d, next_image_x_d, i, dust_params, problem_params);
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

    sprintf(fileName, "%s/cells_gas_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_gas_T = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_gas_T, "%lf %lf %lf\n",prev_x_g[i], prev_gvelocity[i],prev_grho[i]);
    }
    fclose(xy_gas_T);

    sprintf(fileName, "%s/cells_dust_T_d2g=%lf_h=%lf_tau=%lf_T=%lf_K=%lf.dat", DATA_DIR,
            problem_params.d2g, problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_dust_T = fopen(fileName, "w");
    for (int i = 0; i < damount; ++i)
    {
        fprintf(xy_dust_T, "%lf %lf %lf\n", prev_x_d[i], prev_dvelocity[i], prev_drho[i]);
    }
    fclose(xy_dust_T);
}
