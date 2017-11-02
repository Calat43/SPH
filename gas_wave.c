#include "gas_wave.h"

//начальное распределение плотности газа
double gdensity_distribution(double x)
{
    return sin(2*pi*x)/100. + 1;
}

//начальное распределение скорости газа
double gvelocity_distribution(double x)
{
    return sin(2.*pi*x)/100.;
}

//масса газа, находящаяся при постоянной плотности из предположения, что масса всех частиц одинакова
double found_flat_gmass(double * image_x_g, double density, int i, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    int hN = (int)floor(problem_params.h * amount);

    double mass = 0;

    for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
        // for (int j = 0; j < 3 * amount - 2; ++j)
    {
        mass += spline_kernel(image_x_g[amount], image_x_g[j], problem_params);
    }
    return density / mass;
}

//заполнение массива, содержащего массы частиц газа
void fill_gmass(double * gmass, double * x_g, double * image_x_g, double average_grho, int i,
                ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    for(int k = 0; k < amount; ++k)
    {
        //gmass[i] = 1/ (double)amount;
        gmass[k] = gdensity_distribution(x_g[k]) / average_grho * found_flat_gmass(image_x_g, average_grho, i,
                                                                                   particle_params, problem_params);
    }
}

void fill_initial_gvelocity(double * gvelocity, double * x_g, ParticleParams params)
{
    for (int i = 0; i < params.amount; ++i)
    {
        gvelocity[i] = gvelocity_distribution(x_g[i]);
    }
}

double found_next_image_gvelocity(double prev_grho, double * prev_image_grho, double * image_gmass, double * x_g, double * image_x_g,
                                  double prev_gvelocity, int i,  ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double c_s = problem_params.c_s;

    double velocity = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        velocity += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(x_g[i], image_x_g[j], problem_params);
    }
    return -tau * pow(c_s, 2) * velocity + prev_gvelocity;
}

void only_gas_wave(ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    //Блок массивов для газа.BEGIN
    double prev_x_g[amount];
    double next_x_g[amount];

    double gmass[amount];

    double prev_gvelocity[amount];
    double next_gvelocity[amount];
    double prev_grho[amount];
    double next_grho[amount];

    double prev_image_x_g[3*amount - 2];
    double next_image_x_g[3*amount - 2];

    double image_gmass[3*amount - 2];

    double prev_image_gvelocity[3*amount - 2];
    double next_image_gvelocity[3*amount - 2];
    double prev_image_grho[3*amount - 2];
    double next_image_grho[3*amount - 2];

    coordinate_distribution(prev_x_g, particle_params);
    fill_image_x(prev_image_x_g, particle_params);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, 0, particle_params, problem_params);
    fill_image(image_gmass, gmass, particle_params);

    fill_initial_rho(prev_grho, image_gmass, prev_x_g, prev_image_x_g, particle_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, particle_params);

    fill_image(prev_image_grho, prev_grho, particle_params);
    fill_image(prev_image_gvelocity, prev_gvelocity, particle_params);
    //Блок массивов для газа.END

    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = 0; i < amount; ++i)
    {
        fprintf(fout, "%lf %lf\n", prev_x_g[i], prev_grho[i]);
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
            fprintf(rho_frame, "%lf %0.15lf\n", prev_x_g[i], prev_grho[i]);
        }
        fclose(rho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/velocity/frame_%d.dat", frameId);
        FILE * velocity_frame = fopen(fileName, "w");
        for (int i = 0; i < amount; ++i)
        {
            fprintf(velocity_frame, "%lf %0.15lf\n", prev_x_g[i], prev_gvelocity[i]);
        }
        fclose(velocity_frame);

        for(int i = 0; i < amount; ++i)
        {
            next_grho[i] = found_next_rho(image_gmass, prev_x_g, prev_image_x_g, i, particle_params, problem_params);
            next_gvelocity[i] = found_next_image_gvelocity(prev_grho[i], prev_image_grho, image_gmass, prev_x_g,
                                                           prev_image_x_g, prev_gvelocity[i], i, particle_params, problem_params);
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }

        fill_image(next_image_grho, next_grho, particle_params);
        fill_image(next_image_gvelocity, next_gvelocity, particle_params);

        for (int i = 0; i < 3*amount - 2; ++i)
        {
            next_image_x_g[i] = found_next_coordinate(prev_image_x_g[i], prev_image_gvelocity[i], problem_params);
        }

        for(int i = 0; i < amount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];

        }

        for(int i = 0; i < 3 * amount - 2; ++i)
        {
            prev_image_grho[i] = next_image_grho[i];
            prev_image_gvelocity[i] = next_image_gvelocity[i];
            prev_image_x_g[i] = next_image_x_g[i];

        }
    }
}