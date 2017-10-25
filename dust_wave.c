#include "dust_wave.h"

//начальное распределение плотности пыли
double ddensity_distribution(double x)
{
    return sin(2*pi*x)/100 + 0.011;
    //return sin(2*pi*x)/100. + 1;
}

//начальное распределение скорости пыли
double dvelocity_distribution(double x)
{
    return sin(2.*pi*x)/100.;
}

//масса пыли, находящаяся при постоянной плотности из предположения, что масса всех частиц одинакова
double found_flat_dmass(double * image_x_d, double density, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double mass = 0;
    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        mass += spline_kernel(image_x_d[amount], image_x_d[j], problem_params);
    }
    return density / mass;
}

//заполнение массива, содержащего массы частиц пыли
void fill_dmass(double * dmass, double * x_d, double * image_x_d, double average_drho,
                ParticleParams particle_params, ProblemParams problem_params)
{
    double amount = particle_params.amount;

    for(int i = 0; i < amount; ++i)
    {
        //gmass[i] = 1/ (double)amount;
        dmass[i] = ddensity_distribution(x_d[i]) / average_drho * found_flat_dmass(image_x_d, average_drho,
                                                                                   particle_params, problem_params);
    }
}

void fill_initial_dvelocity(double * dvelocity, double * x_d, ParticleParams params)
{
    for (int i = 0; i < params.amount; ++i)
    {
        dvelocity[i] = dvelocity_distribution(x_d[i]);
    }
}

double found_next_image_dvelocity(double prev_dvelocity)
{
    return prev_dvelocity;
}

void only_dust_wave(ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    //Блок массивов для пыли.BEGIN
    double prev_x_d[amount];
    double next_x_d[amount];

    double dmass[amount];

    double prev_dvelocity[amount];
    double next_dvelocity[amount];
    double prev_drho[amount];
    double next_drho[amount];

    double prev_image_x_d[3 * amount - 2];
    double next_image_x_d[3 * amount - 2];

    double image_dmass[3 * amount - 2];

    double prev_image_dvelocity[3 * amount - 2];
    double next_image_dvelocity[3 * amount - 2];
    double prev_image_drho[3 * amount - 2];
    double next_image_drho[3 * amount - 2];

    coordinate_distribution(prev_x_d, particle_params);
    fill_image_x(prev_image_x_d, particle_params);

    double average_drho = ddensity_distribution(0);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, particle_params, problem_params);
    fill_image(image_dmass, dmass, particle_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, particle_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, particle_params);

    fill_image(prev_image_drho, prev_drho, particle_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, particle_params);
    //Блок массивов для пыли.END

    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = 0; i < amount; ++i)
    {
        fprintf(fout, "%lf %lf\n", prev_x_d[i], prev_drho[i]);
    }

    double T = problem_params.T;
    double tau = problem_params.tau;

    char fileName[512];
    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        sprintf(fileName, "/home/calat/CLionProjects/particles/rho/frame_%d.dat", frameId);
        FILE * rho_frame = fopen(fileName, "w");
        for (int i = 0; i < amount; ++i)
        {
            fprintf(rho_frame, "%lf %0.15lf\n", prev_x_d[i], prev_drho[i]);
        }
        fclose(rho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/velocity/frame_%d.dat", frameId);
        FILE * velocity_frame = fopen(fileName, "w");
        for (int i = 0; i < amount; ++i)
        {
            fprintf(velocity_frame, "%lf %0.15lf\n", prev_x_d[i], prev_dvelocity[i]);
        }
        fclose(velocity_frame);

        for(int i = 0; i < amount; ++i)
        {
            next_drho[i] = found_next_rho(image_dmass, prev_x_d, prev_image_x_d, i, particle_params, problem_params);
            next_dvelocity[i] = found_next_image_dvelocity(prev_dvelocity[i]);
            next_x_d[i] = found_next_coordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }

        fill_image(next_image_drho, next_drho, particle_params);
        fill_image(next_image_dvelocity, next_dvelocity, particle_params);

        for (int i = 0; i < 3 * amount - 2; ++i)
        {
            next_image_x_d[i] = found_next_coordinate(prev_image_x_d[i], prev_image_dvelocity[i], problem_params);
        }

        for(int i = 0; i < amount; ++i)
        {
            prev_drho[i] = next_drho[i];
            prev_dvelocity[i] = next_dvelocity[i];
            prev_x_d[i] = next_x_d[i];

        }

        for(int i = 0; i < 3 * amount - 2; ++i)
        {
            prev_image_drho[i] = next_image_drho[i];
            prev_image_dvelocity[i] = next_image_dvelocity[i];
            prev_image_x_d[i] = next_image_x_d[i];

        }
    }
}