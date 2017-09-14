#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

const double pi = 3.14159265358;

const double T = 0.5;
const double left = 0;
const double right = 1;
const int amount = 400;
const double h = 0.01;
const double tau = 0.001;
const double c_s = 1;

//распределение плотности частиц в пространстве
double particle_distribution(double x)
{
    return 1;
            //sin(2.*pi*x)/100. + 1;
}

//начальное распределение скорости
double velocity_distribution(double x)
{
    return 0;
            //sin(2.*pi*x)/100.;
}

double gauss_kernel(double x_a, double x_b)
{
    return 1 / h / sqrt(pi) * exp(-pow(x_a - x_b, 2) / h / h);
}

double gauss_gradient(double x_a, double x_b)
{
    return - 2 * (x_a - x_b) / pow(h, 2) * gauss_kernel(x_a, x_b);
}

double spline_kernel(double x_a, double x_b)
{
    double r = fabs(x_a - x_b);
    double result = 0;
    if (r / h >= 0 && r / h <= 1)
    {
        result = 1 - 3. / 2 * pow(r, 2) + 3. / 4 * pow(r, 3);
        return 1./h/pi * result;
    }
    if (r / h >= 1 && r / h <= 2)
    {
        result = 1. / 4 * pow((2. - r), 3);
        return 1./h/pi * result;
    }
    return 1./h/pi * result;
}

double spline_gradient0(double x_a, double x_b)
{
    double r = x_a - x_b;
    double abs_r = fabs(x_a - x_b);
    double result = 0;
    if (abs_r / h >= 0 && abs_r / h <= 1)
    {
        result = - 3 * r + 9. / 4 * pow(r, 2);
        return 1./h/pi * result;
    }
    if (abs_r / h >= 1 && abs_r / h <= 2)
    {
        result = - 3. / 4 * pow((2. - r), 2);
        return 1./h/pi * result;
    }
    return 1./h/pi * result;
}

double derivative_r(double x_a, double x_b)
{
    if (x_a > x_b)
    {
        return 1;
    }
    else if (x_a < x_b)
    {
        return -1;
    }
    return 0;
}

double spline_gradient(double x_a, double x_b)
{
    double r = fabs(x_a - x_b);
    double result = 0;
    if (r / h >= 0 && r / h <= 1)
    {
        if (x_a >= x_b)
        {
            result = -3 * (x_a - x_b) + 9. / 4 * pow(x_a - x_b, 2);
        }
        else
        {
            result = -3 * (x_a - x_b) - 9. / 4 * pow(x_a - x_b, 2);
        }
    }
    else if (r / h >= 1 && r / h <= 2)
    {
        if (x_a >= x_b)
        {
            result = - 3. / 4 * pow((2. - (x_a - x_b)), 2);
        }
        else
        {
            result = 3. / 4 * pow((2. + (x_a - x_b)), 2);
        }
    }
    return 1./h/pi * result;
}

//координаты точек располагаются равномерно (пока что)
void fill_x(double * x)
{
    for (int i = 0; i < amount; ++i)
    {
        x[i] = left + i * (right - left) / (amount - 1);
    }
}

double fill_image_x(double * image_x)
{
    double new_left = left - (right - left);
    double new_right = right + (right + left);

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        image_x[j] = new_left + j * (new_right - new_left) / (3*(amount - 1));
    }
}

void fill_mass(double * mass, double * x, double average, double flat_mass)
{
    for(int i = 0; i < amount; ++i)
    {
        //mass[i] = particle_distribution(x[i]);
        mass[i] = particle_distribution(x[i]) / average * flat_mass;
    }
}

void fill_image_mass(double * image_mass, double * image_x, double avarage_x, double flat_mass)
{
    for (int i = 0; i < amount; ++i)
    {
        double real_mass_i = 1./(double)amount;
                //particle_distribution(image_x[amount - 1 + i]) / avarage_x * flat_mass;
        image_mass[i] = real_mass_i;

        image_mass[amount - 1 + i] = 0;
        image_mass[2 * amount - 2 + i] = 0;

        //image_mass[amount - 1 + i] = real_mass_i;
        //image_mass[2 * amount - 2 + i] = real_mass_i;
    }
}

double found_flat_mass(double * image_x, double density)
{
    double mass = 0;
    {
        for (int j = 0; j < 3*amount - 2; ++j)
        {
            mass += spline_kernel(image_x[amount - 1], image_x[j]);
        }
    }
    return density / mass;
}

void fill_flat_rho(double * rho, double * x, double * image_x, double flat_mass)
{
    for (int i = 0; i < amount; ++i)
    {
        rho[i] = 0;
        for (int j = 0; j < 3 * amount - 2; ++j)
        {
            rho[i] += flat_mass * spline_kernel(x[i], image_x[j]);
        }
    }
}

void fill_rho(double * image_rho, double  * image_mass, double * image_x)
{
    for(int i = 0; i < amount; ++i)
    {
        double real_rho_i = 0;
        for (int j = 0; j < 3*amount - 2; ++j)
        {
            real_rho_i += image_mass[j] * spline_kernel(image_x[amount - 1 + i], image_x[j]);
        }
        image_rho[i] = real_rho_i;
        image_rho[amount - 1 + i] = real_rho_i;
        image_rho[2*amount - 2 + i] = real_rho_i;
    }
}

void fill_image_rho(double * image_rho, double  * image_mass, double * image_x)
{
    for(int i = 0; i < amount; ++i)
    {
        double real_rho_i = 0;
        for (int j = 0; j < 3*amount - 2; ++j)
        {
            real_rho_i += image_mass[j] * spline_kernel(image_x[amount - 1 + i], image_x[j]);
        }
        image_rho[i] = real_rho_i;
        image_rho[amount - 1 + i] = real_rho_i;
        image_rho[2*amount - 2 + i] = real_rho_i;
    }
}

double fill_velocity(double * velocity, double * x)
{
    for (int i = 0; i < amount; ++i)
    {
        velocity[i] = velocity_distribution(x[i]);
    }
}

double fill_image_velocity(double * image_velocity, double * image_x)
{
    for(int i = 0; i < amount; ++i)
    {
        double real_v_i = velocity_distribution(image_x[amount - 1 + i]);
        image_velocity[i] = real_v_i;
        image_velocity[amount - 1 + i] = real_v_i;
        image_velocity[2*amount - 2 + i] = real_v_i;
    }
}

// TODO replace with imaginary version
double rho_error(double * rho, double * x)
{
    double error = 0;
    for (int i = 0; i < amount; ++i)
    {
        if(fabs(rho[i] - particle_distribution(x[i])) > error)
        {
            error = fabs(rho[i] - particle_distribution(x[i]));
        }
    }
    return error;
}


double found_next_rho(double * image_velocity, double * image_mass, double * image_x, double prev_rho, double tau, int i)
{
    double rho = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        rho += image_mass[j] * (image_velocity[i] - image_velocity[j]) * spline_gradient(image_x[i], image_x[j]);
    }
    return tau * rho + prev_rho;
}

double found_next_gvelocity(double * image_rho, double * image_mass, double * image_x, double prev_velocity, double c_s, double tau, int i)
{
    double velocity = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        velocity += image_mass[j] * (1. / image_rho[j] + 1. / image_rho[i]) * spline_gradient(image_x[i], image_x[j]);
    }
    return -tau * pow(c_s, 2) * velocity + prev_velocity;
}

int main()
{
    /*   left image   real   right image
     * |............|......|.............|
     * ^            ^      ^             ^
     * 0    (amount-1)    (2*amount-2)  (3*amount-3)
     */
    double image_x[3 * amount - 2];
    double image_rho[3 * amount - 2];
    double image_mass[3 * amount - 2];

    double prev_gvelocity[3*amount - 2];
    double next_gvelocity[3*amount - 2];
    double prev_grho[3*amount - 2];
    double next_grho[3*amount - 2];

    double x[amount];
    double mass[amount];

    fill_x(x);

    //fill_image_x(image_x);
    //fill_image_velocity(prev_gvelocity, image_x);

    double average_f = particle_distribution(0);
    double flat_mass = found_flat_mass(image_x, average_f);

    fill_mass(mass, x, average_f, flat_mass);
    fill_image_mass(image_mass, image_x, average_f, flat_mass);

    fill_image_rho(image_rho, image_mass, image_x);

    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = amount - 1; i < 2 * amount - 1; ++i)
    {
        fprintf(fout, "%lf %lf\n", image_x[i], image_mass[i]);
    }

    for(int i = 0; i < 3 * amount - 2; ++i)
    {
        prev_grho[i] = image_rho[i];
    }

    char fileName[512];
    for (int frameId = 0; frameId < floor(T / tau); ++frameId)
    {
        sprintf(fileName, "/home/calat/CLionProjects/particles/rho/frame_%d.dat", frameId);
        FILE * rho_frame = fopen(fileName, "w");
        for (int i = 0; i < 3 * amount - 2; ++i)
        {
            fprintf(rho_frame, "%lf %0.15lf\n", image_x[i], prev_grho[i]);
        }
        fclose(rho_frame);

        sprintf(fileName, "/home/calat/CLionProjects/particles/velocity/frame_%d.dat", frameId);
        FILE * velocity_frame = fopen(fileName, "w");
        for (int i = 0; i < 3 * amount - 2; ++i)
        {
            fprintf(velocity_frame, "%lf %0.15lf\n", image_x[i], prev_gvelocity[i]);
        }
        fclose(velocity_frame);

        for(int i = 0; i < 3 * amount - 2; ++i)
        {
            next_grho[i] = found_next_rho(prev_gvelocity, image_mass, image_x, prev_grho[i], tau, i);
            next_gvelocity[i] = found_next_gvelocity(prev_grho, image_mass, image_x, prev_gvelocity[i], c_s, tau, i);
        }

        for(int i = 0; i < 3 * amount - 2; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
        }
    }

    return 0;
}