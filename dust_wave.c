#include "dust_wave.h"

//начальное распределение плотности пыли
double ddensity_distribution(double x)
{
    return sin(2*pi*x)/100. + 1;
}

//начальное распределение скорости пыли
double dvelocity_distribution(double x)
{
    return sin(2.*pi*x)/100.;
}

//начальное положение частиц пыли (заполняется сразу как массив, равномерно)
void dcoordinate_distribution(double * x_d, ParticleParams params)
{
    int amount = params.amount;
    double left = params.left;
    double right = params.right;

    for (int i = 0; i < amount; ++i)
    {
        x_d[i] = left + i * (right - left) / (amount - 1);
    }
}

//начальное положение частиц пыли (в том числе мнимых)
double fill_dimage_x(double * image_x, ParticleParams params)
{
    int amount = params.amount;
    double left = params.left;
    double right = params.right;

    double new_left = left - (right - left);
    double new_right = right + (right + left);

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        image_x[j] = new_left + j * (new_right - new_left) / (3*(amount - 1));
    }
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

void fill_dimage(double * image, double * real, ParticleParams params)
{
    int amount = params.amount;

    for (int i = 0; i < amount; ++i)
    {
        image[i] = real[i];
        image[amount - 1 + i] = real[i];
        image[2 * amount - 2 + i] = real[i];
    }
}


//заполнение массива начального распределения плотности пыли
void fill_initial_drho(double * drho, double  * image_dmass, double * x_d, double * image_x_d,
                       ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    for(int i = 0; i < amount; ++i)
    {
        for (int j = 0; j < 3 * amount - 2; ++j)
        {
            drho[i] += image_dmass[j] * spline_kernel(x_d[i], image_x_d[j], problem_params);
        }
    }
}

double fill_initial_dvelocity(double * dvelocity, double * x_d, ParticleParams params)
{
    for (int i = 0; i < params.amount; ++i)
    {
        dvelocity[i] = dvelocity_distribution(x_d[i]);
    }
}

double found_next_drho(double * dmass, double * x_d, int i, ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double rho = 0;
    for(int j = 0; j < amount; ++j)
    {
        rho += dmass[j] * spline_kernel(x_d[i], x_d[j], problem_params);
    }
    return rho;
}

double found_next_image_drho(double * image_dmass, double * x_d, double * image_x_d, int i,
                             ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double rho = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        rho += image_dmass[j] * spline_kernel(x_d[i], image_x_d[j], problem_params);
    }
    return rho;
}

double found_next_dvelocity(double prev_dvelocity)
{
    return prev_dvelocity;
}

double found_next_image_dvelocity(double prev_dvelocity)
{
    return prev_dvelocity;
}

double found_next_dcoordinate(double prev_x_d, double prev_dvelocity, ProblemParams params)
{
    return prev_x_d + params.tau * prev_dvelocity;
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

    dcoordinate_distribution(prev_x_d, particle_params);
    fill_dimage_x(prev_image_x_d, particle_params);

    double average_drho = ddensity_distribution(0);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, particle_params, problem_params);
    fill_dimage(image_dmass, dmass, particle_params);

    fill_initial_drho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, particle_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, particle_params);

    fill_dimage(prev_image_drho, prev_drho, particle_params);
    fill_dimage(prev_image_dvelocity, prev_dvelocity, particle_params);
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
            next_drho[i] = found_next_image_drho(image_dmass, prev_x_d, prev_image_x_d, i, particle_params, problem_params);
            next_dvelocity[i] = found_next_image_dvelocity(prev_dvelocity[i]);
            next_x_d[i] = found_next_dcoordinate(prev_x_d[i], prev_dvelocity[i], problem_params);
        }

        fill_dimage(next_image_drho, next_drho, particle_params);
        fill_dimage(next_image_dvelocity, next_dvelocity, particle_params);

        for (int i = 0; i < 3 * amount - 2; ++i)
        {
            next_image_x_d[i] = found_next_dcoordinate(prev_image_x_d[i], prev_image_dvelocity[i], problem_params);
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