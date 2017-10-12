#include "common_use.h"

const double pi = 3.14159265358;

double spline_kernel(double x_a, double x_b, ProblemParams params)
{
    double h = params.h;
    double r = fabs(x_a - x_b);
    double q = r / h;
    double result = 0;
    if (r / h >= 0 && r / h <= 1)
    {
        result = 1 - 3. / 2 * pow(q, 2) + 3. / 4 * pow(q, 3);
        return 2./ 3. / h * result;
    }
    if (r / h >= 1 && r / h <= 2)
    {
        result = 1. / 4 * pow((2. - q), 3);
        return 2./ 3. / h * result;
    }
    return 2./ 3. / h * result;
}

double spline_gradient(double x_a, double x_b, ProblemParams params)
{
    double h = params.h;
    double r = fabs(x_a - x_b);
    double q = r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }
    if (x_a > x_b)
    {
        return 2./3. / h / h * result;
    }
    if (x_a == x_b)
    {
        return 0;
    }
    if (x_a < x_b)
    {
        return - 2./3. / h / h * result;
    }
}

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params)
{
    return prev_x + params.tau * prev_vel;
}

void coordinate_distribution(double * x, ParticleParams params)
{
    int amount = params.amount;
    double left = params.left;
    double right = params.right;

    for (int i = 0; i < amount; ++i)
    {
        x[i] = left + i * (right - left) / (amount - 1);
    }
}

//начальное положение частиц (в том числе мнимых)
double fill_image_x(double * image_x, ParticleParams params)
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

void fill_image(double * image, double * real, ParticleParams params)
{
    int amount = params.amount;

    for (int i = 0; i < amount; ++i)
    {
        image[i] = real[i];
        image[amount - 1 + i] = real[i];
        image[2 * amount - 2 + i] = real[i];
    }
}

//заполнение массива начального распределения плотности
void fill_initial_rho(double * rho, double  * image_mass, double * x, double * image_x,
                       ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    for(int i = 0; i < amount; ++i)
    {
        for (int j = 0; j < 3 * amount - 2; ++j)
        {
            rho[i] += image_mass[j] * spline_kernel(x[i], image_x[j], problem_params);
        }
    }
}

double found_next_rho(double * image_mass, double * x, double * image_x, int i,
                             ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;

    double rho = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        rho += image_mass[j] * spline_kernel(x[i], image_x[j], problem_params);
    }
    return rho;
}

double point_value(double x, double * image_function, double * image_mass, double * image_rho, double * image_x,
                   ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double result = 0;

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_mass[j] * image_function[j] / image_rho[j] * spline_kernel(x, image_x[j], problem_params);
    }

    return result;
}