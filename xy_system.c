#include "xy_system.h"

//заполнение массива x, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void compute_x(double * x, double * gvelocity, double * dvelocity, ParticleParams params)
{
    for(int i = 0; i < params.amount; ++i)
    {
        x[i] = gvelocity[i] - dvelocity[i];
    }
}

//заполнение массива y, опирающееся на предположение о том, что начальные данные заданы в одних и тех же точках
void compute_y(double * y, double * gvelocity, double * dvelocity, double * grho, double * drho, ParticleParams params)
{
    for(int i = 0; i < params.amount; ++i)
    {
        y[i] = gvelocity[i] + drho[i] / grho[i] * dvelocity[i];
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
    int hN = (int)floor(problem_params.h * amount);

    double result = 0;
    double sum = 0;
    double denom = 0;
    //for(int j = 0; j < 3 * amount - 2; ++j)
    for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(r, prev_image_x_g[j],
                                                                                             problem_params);
    }
    denom = 1. + tau * K * (1. / prev_grho + 1. / prev_drho);
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

    //for(int j = 0; j < 3 * amount - 2; ++j)
    for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    {
        sum += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(r, prev_image_x_g[j],
                                                                                             problem_params);
    }
    result = - tau * pow(c_s, 2) * sum + prev_y;

    return result;
}

double found_next_dvelocity(double next_x, double next_y, double rho_div)
{
    return (next_y - next_x) / (1. + rho_div);
}

double found_next_gvelocity(double next_x, double next_y, double rho_div)
{
    return (next_y + rho_div * next_x)  / (1 + rho_div);
}

double near_velocity(double r, double * coordinate, double * velocity, ParticleParams params)
{
    int nearly = 0;

    for (int i = 0; i < params.amount; ++i)
    {
        if (fabs(r - coordinate[i]) < fabs(r - coordinate[nearly]))
        {
            nearly = i;
        }
    }

    return velocity[nearly];
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
    double next_image_grho[3 * gamount - 2];

    coordinate_distribution(prev_x_g, gas_params);
    fill_image_x(prev_image_x_g, gas_params);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, 0, gas_params, problem_params);
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
    double next_image_drho[3 * damount - 2];

    coordinate_distribution(prev_x_d, dust_params);
    fill_image_x(prev_image_x_d, dust_params);

    double average_drho = ddensity_distribution(0);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, 0, dust_params, problem_params);
    fill_image(image_dmass, dmass, dust_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params);

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

    for(int i = 0; i < gamount; ++i)
    {
        prev_grho_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_grho, dust_params);
        prev_drho_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_drho, gas_params);
        //prev_drho_xg[i] = point_value_for_rho(prev_x_g[i], image_dmass, prev_image_x_d, dust_params, problem_params);
        //prev_grho_xd[i] = point_value_for_rho(prev_x_d[i], image_gmass, prev_image_x_g, gas_params, problem_params);

        prev_gvel_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_gvelocity, dust_params);
        prev_dvel_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_dvelocity, gas_params);
        //prev_gvel_xd[i] = point_value(prev_x_d[i], prev_image_gvelocity, image_gmass, prev_image_grho, prev_image_x_g, gas_params, problem_params);
        //prev_dvel_xg[i] = point_value(prev_x_g[i], prev_image_dvelocity, image_dmass, prev_image_drho, prev_image_x_d, dust_params, problem_params);
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

    compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
    compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params);

    compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
    compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params);
    //Массивы для x,y.END


    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = 0; i < gamount; ++i)
    {
        fprintf(fout, "%lf\n", prev_dvelocity[i] - prev_dvel_xg[i]);
    }


    char fileName[512];

    sprintf(fileName, "/home/calat/CLionProjects/particles/xy_h=%lf_tau=%lf_T=%lf_K=%lf.dat", problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_0 = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_0, "%lf %lf %lf %lf %lf\n", prev_x_g[i], prev_dvelocity[i], prev_gvelocity[i], prev_drho[i], prev_grho[i]);
    }
    fclose(xy_0);

    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        printf("%d\n", frameId);
    /*
        sprintf(fileName, "/home/calat/CLionProjects/particles/delta_vel/delta_gvel_%d.dat", frameId);
        FILE * delta_gvel_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(delta_gvel_frame, "%lf %lf %lf %lf \n", prev_x_g[i], prev_gvelocity[i], prev_x_d[i], prev_gvel_xd[i]);
        }
        fclose(delta_gvel_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/delta_vel/delta_dvel_%d.dat", frameId);
        FILE * delta_dvel_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(delta_dvel_frame, "%lf %lf %lf %lf \n", prev_x_d[i], prev_dvelocity[i], prev_x_g[i], prev_dvel_xg[i]);
        }
        fclose(delta_dvel_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/delta_vel/delta_grho_%d.dat", frameId);
        FILE * delta_grho_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(delta_grho_frame, "%lf %lf %lf %lf \n", prev_x_g[i], prev_grho[i], prev_x_d[i], prev_grho_xd[i]);
        }
        fclose(delta_grho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/delta_vel/delta_drho_%d.dat", frameId);
        FILE * delta_drho_frame = fopen(fileName, "w");
        for (int i = 0; i < gamount; ++i)
        {
            fprintf(delta_drho_frame, "%lf %lf %lf %lf \n", prev_x_d[i], prev_drho[i], prev_x_g[i], prev_drho_xg[i]);
        }
        fclose(delta_drho_frame);

     */
/*
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

        sprintf(fileName, "/home/calat/CLionProjects/particles/xgas/xgas_%d.dat", frameId);
        FILE * xgas_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(xgas_frame, "%lf %0.15lf\n", prev_x_g[i], prev_x_gas[i]);
        }
        fclose(xgas_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/ygas/ygas_%d.dat", frameId);
        FILE * ygas_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(ygas_frame, "%lf %0.15lf\n", prev_x_g[i], prev_y_gas[i]);
        }
        fclose(ygas_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/xdust/xdust_%d.dat", frameId);
        FILE * xdust_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(xdust_frame, "%lf %0.15lf\n", prev_x_d[i], prev_x_dust[i]);
        }
        fclose(xdust_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/ydust/ydust_%d.dat", frameId);
        FILE * ydust_frame = fopen(fileName, "w");
        for (int i = 0; i < damount; ++i)
        {
            fprintf(ydust_frame, "%lf %0.15lf\n", prev_x_d[i], prev_y_dust[i]);
        }
        fclose(ydust_frame);
*/


        for(int i = 0; i < gamount; ++i)
        {
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }
        for(int j = 0; j < 3 * gamount - 2; ++j)
        {
            next_image_x_d[j] = found_next_coordinate(prev_image_x_d[j], prev_image_dvelocity[j], problem_params);
            next_image_x_g[j] = found_next_coordinate(prev_image_x_g[j], prev_image_gvelocity[j], problem_params);
        }
        for(int i = 0; i < gamount; ++i)
        {
            next_drho[i] = found_next_rho(image_dmass, next_x_d, next_image_x_d, i, dust_params, problem_params);
            next_grho[i] = found_next_rho(image_gmass, next_x_g, next_image_x_g, i, dust_params, problem_params);
        }

        for(int i = 0; i < gamount; ++i)
        {
            next_grho_xd[i] = near_velocity(next_x_d[i], next_x_g, next_grho, dust_params);
            next_drho_xg[i] = near_velocity(next_x_g[i], next_x_d, next_drho, gas_params);
            //next_drho_xg[i] = point_value_for_rho(next_x_g[i], image_dmass, next_image_x_d, dust_params, problem_params);
            //next_grho_xd[i] = point_value_for_rho(next_x_d[i], image_gmass, next_image_x_g, gas_params, problem_params);
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
                                     prev_drho[i], prev_x_d[i], i, dust_params, problem_params);
            next_y_dust[i] = found_next_y(prev_y_dust[i], image_gmass, prev_image_grho, prev_image_x_g, prev_grho_xd[i],
                                     prev_x_d[i], i, dust_params, problem_params);
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

            prev_drho[i] = next_drho[i];
            prev_dvelocity[i] = next_dvelocity[i];
            prev_x_d[i] = next_x_d[i];

            prev_x_gas[i] = next_x_gas[i];
            prev_y_gas[i] = next_y_gas[i];

            prev_x_dust[i] = next_x_dust[i];
            prev_y_dust[i] = next_y_dust[i];

            prev_drho_xg[i] = next_drho_xg[i];
            prev_grho_xd[i] = next_grho_xd[i];
        }

        for(int i = 0; i < 3 * gamount - 2; ++i)
        {
            prev_image_grho[i] = next_image_grho[i];
            prev_image_gvelocity[i] = next_image_gvelocity[i];
            prev_image_x_g[i] = next_image_x_g[i];

            prev_image_drho[i] = next_image_drho[i];
            prev_image_dvelocity[i] = next_image_dvelocity[i];
            prev_image_x_d[i] = next_image_x_d[i];

        }

        for(int i = 0; i < gamount; ++i)
        {
            prev_gvel_xd[i] = near_velocity(prev_x_d[i], prev_x_g, prev_gvelocity, dust_params);
            prev_dvel_xg[i] = near_velocity(prev_x_g[i], prev_x_d, prev_dvelocity, gas_params);
            //prev_gvel_xd[i] = point_value(prev_x_d[i], prev_image_gvelocity, image_gmass, prev_image_grho, prev_image_x_g, gas_params, problem_params);
            //prev_dvel_xg[i] = point_value(prev_x_g[i], prev_image_dvelocity, image_dmass, prev_image_drho, prev_image_x_d, dust_params, problem_params);
        }

        compute_x(prev_x_gas, prev_gvelocity, prev_dvel_xg, gas_params);
        compute_y(prev_y_gas, prev_gvelocity, prev_dvel_xg, prev_grho, prev_drho_xg, gas_params);

        compute_x(prev_x_dust, prev_gvel_xd, prev_dvelocity, dust_params);
        compute_y(prev_y_dust, prev_gvel_xd, prev_dvelocity, prev_grho_xd, prev_drho, dust_params);
    }

    sprintf(fileName, "/home/calat/CLionProjects/particles/xy_T_h=%lf_tau=%lf_T=%lf_K=%lf.dat", problem_params.h, problem_params.tau, problem_params.T, problem_params.K);
    FILE * xy_T = fopen(fileName, "w");
    for (int i = 0; i < gamount; ++i)
    {
        fprintf(xy_T, "%lf %lf %lf %lf %lf\n", prev_x_g[i], prev_dvelocity[i], prev_gvelocity[i], prev_drho[i], prev_grho[i]);
    }
    fclose(xy_T);

}