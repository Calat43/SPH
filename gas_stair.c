#include "gas_stair.h"
#include "common_use.h"

//начальное распределение плотности газа
double stair_gdensity_distribution(double x)
{
    return 1;
            //sin(2*pi*x)/100. + 1;
}

//начальное распределение скорости газа
double stair_gvelocity_distribution(double x)
{
    return 0;
    //sin(2.*pi*x)/100.;
}

//начальное положение частиц газа (заполняется сразу как массив, равномерно)
void stair_gcoordinate_distribution(double * x_g, ParticleParams params)
{
    int amount = params.amount;
    double left = params.left;
    double right = params.right;

    for (int i = 0; i < amount; ++i)
    {
        x_g[i] = left + i * (right - left) / (amount - 1);
    }
}

//масса газа, находящаяся при постоянной плотности из предположения, что масса всех частиц одинакова
double stair_found_flat_gmass(double * x_g, double density, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double mass = 0;
    for (int j = 0; j < amount; ++j)
    {
        mass += spline_kernel(x_g[amount], x_g[j], problem_params);
    }
    return density / mass;
}

//заполнение массива, содержащего массы частиц газа
void stair_fill_gmass(double * gmass, double * x_g, double average_grho, ParticleParams params)
{
    int amount = params.amount;

    for(int i = 0; i < amount; ++i)
    {
        gmass[i] = 1/ (double)amount;
        //gmass[i] = gdensity_distribution(x_g[i]) / average_grho * found_flat_gmass(x_g, average_grho);
    }
}

//заполнение массива начального распределения плотности газа
void stair_fill_initial_grho(double * grho, double  * gmass, double * x_g, ParticleParams particle_params,
                       ProblemParams problem_params)
{
    int amount = particle_params.amount;

    for(int i = 0; i < amount; ++i)
    {
        for (int j = 0; j < amount; ++j)
        {
            grho[i] += gmass[j] * spline_kernel(x_g[i], x_g[j], problem_params);
        }
    }
}

double stair_fill_initial_velocity(double * gvelocity, double * x_g, ParticleParams params)
{

    for (int i = 0; i < params.amount; ++i)
    {
        gvelocity[i] = stair_gvelocity_distribution(x_g[i]);
    }
}

double stair_found_next_grho(double * gmass, double * x_g, int i, ParticleParams particle_params,
                       ProblemParams problem_params)
{
    double rho = 0;
    for(int j = 0; j < particle_params.amount; ++j)
    {
        rho += gmass[j] * spline_kernel(x_g[i], x_g[j], problem_params);
    }
    return rho;
}

double stair_found_next_gvelocity(double * prev_grho, double * gmass, double * x_g, double prev_gvelocity, int i,
                            ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double tau = problem_params.tau;
    double c_s = problem_params.c_s;

    double velocity = 0;
    for(int j = 0; j < amount; ++j)
    {
        velocity += gmass[j] * (1. / prev_grho[j] + 1. / prev_grho[i]) * spline_gradient(x_g[i], x_g[j], problem_params);
    }
    return -tau * pow(c_s, 2) * velocity + prev_gvelocity;
}

double stair_gas_print(ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double prev_x_g[amount];
    double next_x_g[amount];

    double gmass[amount];

    double prev_gvelocity[amount];
    double next_gvelocity[amount];

    double prev_grho[amount];
    double next_grho[amount];

    stair_gcoordinate_distribution(prev_x_g, particle_params);

    double average_grho = stair_gdensity_distribution(0);
    stair_fill_gmass(gmass, prev_x_g, average_grho, particle_params);

    stair_fill_initial_grho(prev_grho, gmass, prev_x_g, particle_params, problem_params);
    stair_fill_initial_velocity(prev_gvelocity, prev_x_g, particle_params);

    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = 0; i < amount; ++i)
    {
        fprintf(fout, "%lf %lf\n", prev_x_g[i], prev_grho[i]);
    }

    char fileName[512];
    for (int frameId = 0; frameId < floor(problem_params.T / problem_params.tau); ++frameId)
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
            next_grho[i] = stair_found_next_grho(gmass, prev_x_g, i, particle_params, problem_params);
            next_gvelocity[i] = stair_found_next_gvelocity(prev_grho, gmass, prev_x_g, prev_gvelocity[i], i,
                                                           particle_params, problem_params);
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i], problem_params);
        }

        for(int i = 0; i < amount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];

        }
    }
}