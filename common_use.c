#include <assert.h>
#include <stdbool.h>
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
    if (r / h > 1 && r / h <= 2)
    {
        result = 1. / 4 * pow((2. - q), 3);
        return 2./ 3. / h * result;
    }
    if(r / h > 2)
    {
        return 0;
    }
    assert(false);
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

    assert(false);
}

double found_next_coordinate(double prev_x, double prev_vel, ProblemParams params)
{
    return prev_x + params.tau * prev_vel;
}
/*
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
void fill_image_x(double * image_x, ParticleParams params)
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
*/
double coord_function(double x_param, double x, int amount, ProblemParams params, ParticleParams particleParams)
{
    double middle = 0;
    if(particleParams.isGas == true)
    {
        middle = params.middle_rho_gas;
    }
    else
    {
        middle = params.d2g;
    }
    return middle * x - params.delta / 2. / pi
                                             * (cos(2*pi*(x_param + x)) - cos(2*pi*x_param)) - middle/amount;
}

double deriv_func(double x_param, double x, ProblemParams params, ParticleParams particleParams)
{
    double middle = 0;
    if(particleParams.isGas == true)
    {
        middle = params.middle_rho_gas;
    }
    else
    {
        middle = params.d2g;
    }
    return middle + params.delta * sin(2*pi*(x_param + x));
}

//ищем расстояние между двумя точками
double newton(double x_param, double x, double exact, int amount, ProblemParams params, ParticleParams particleParams)
{
    double prev_x = x;
    double next_x = 0;
    double tmp = 0;
    int i = 0;

    while(true)
    {
        printf("%d   %0.20lf %0.20lf\n", i, fabs(next_x - prev_x), next_x);
        //tmp = deriv_func(x_param, prev_x, params);
        next_x = prev_x - coord_function(x_param, prev_x, amount, params, particleParams)
                          / deriv_func(x_param, prev_x, params, particleParams);
        if(fabs(next_x - prev_x) < exact)
        {
            return next_x;
        }
        prev_x = next_x;
        i++;
    }
}


void fill_x(double * coord, ProblemParams problemParams, ParticleParams particleParams)
{
    int amount = particleParams.amount;
    double x_estimate = 1. / amount;
    double exact = 1. / 100000000. / amount;
    double prev_x = 0;
    double next_x = 0;

    prev_x = particleParams.left;
    for(int i = 0; i < amount; ++i)
    {
        next_x = prev_x + newton(prev_x,  x_estimate, exact, amount, problemParams, particleParams);
        coord[i] = prev_x + 1./ 2 * newton(prev_x, x_estimate, exact, amount, problemParams, particleParams);
        prev_x = next_x;
    }
}


void fill_image_x(double * image_coord, double * coord, ParticleParams params)
{
    int amount = params.amount;

    for(int i = 0; i < amount; ++i)
    {
        image_coord[i] = coord[i] - 1.;
        image_coord[amount + i] = coord[i];
        image_coord[2 * amount + i] = coord[i] + 1;
    }
}


void dust_fill_image_x(double * image_x, ParticleParams params)
{
    int amount = params.amount;
    double left = params.left;
    double right = params.right;

    double new_left = left - (right - left);
    double new_right = right + (right + left);

    double mix = (right - left) / (amount - 1) / 2.;

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        image_x[j] = new_left + j * (new_right - new_left) / (3*(amount - 1)) + mix;
    }
}

void fill_image(double * image, double * real, ParticleParams params)
{
    int amount = params.amount;

    for (int i = 0; i < amount; ++i)
    {
        assert(!isnan(real[i]));
        image[i] = real[i];
        image[amount + i] = real[i];
        image[2 * amount + i] = real[i];
    }
}

//заполнение массива начального распределения плотности
void fill_initial_rho(double * rho, double  * image_mass, double * x, double * image_x,
                       ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    int hN = (int)floor(problem_params.h * amount);

    for(int i = 0; i < amount; ++i)
    {
        rho[i] = 0;
        //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
        for(int j = 0; j < 3 * amount - 2; ++j)
        {
            rho[i] += image_mass[j] * spline_kernel(x[i], image_x[j], problem_params);
            if(rho == 0)
            {
                assert(false);
            }
        }
    }
}

void fill_mass(double * mass, ParticleParams params, ProblemParams problemParams)
{
    double middle_rho = 0;
    if(params.isGas == true)
    {
        middle_rho = problemParams.middle_rho_gas;
    }
    else
    {
        middle_rho  = problemParams.d2g;
    }

    for(int i = 0; i < params.amount; ++i)
    {
        mass[i] = middle_rho / params.amount;
    }
}

double found_next_rho(double * image_mass, double * x, double * image_x, int i,
                             ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    int hN = (int)floor(problem_params.h * amount);

    double rho = 0;
    //for (int j = i + amount - 2*hN - 1; j < i + amount + 2*hN + 1; ++j)
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        rho += image_mass[j] * spline_kernel(x[i], image_x[j], problem_params);
    }
    return rho;
}

double interpolation_value(double x, double * image_function, double * image_mass, double * image_rho, double * image_x,
                   ParticleParams particle_params, ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double result = 0;
    int hN = (int)floor(problem_params.h * amount);

    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_mass[j] * (image_function[j] / image_rho[j]) * spline_kernel(x, image_x[j], problem_params);
        assert(!isnan(result));
    }

    return result;
}

double image_interpolation_value(double * image_x, double * image_function, double * image_mass, double * image_rho,
                                 double * interpol_x, int i, ParticleParams particleParams, ProblemParams problemParams)
{
    int amount = particleParams.amount;
    double result = 0;
    int hN = (int)floor(problemParams.h * amount);

    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_mass[j] * (image_function[j] / image_rho[j]) * spline_kernel(image_x[i], image_x[j], problemParams);
        assert(!isnan(result));
    }

    return result;
}

double interpolation_value_for_rho(double x, double * image_mass, double * image_x, int i, ParticleParams particle_params,
                           ProblemParams problem_params)
{
    int amount = particle_params.amount;
    double result = 0;
    int hN = (int)floor(problem_params.h * amount);

    //for (int j = i + amount - 1 - 2*hN; j < i + amount + 2*hN + 1; ++j)
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        result += image_mass[j] * spline_kernel(x, image_x[j], problem_params);
    }

    return result;
}