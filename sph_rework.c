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

//начальное положение частиц газа (заполняется сразу как массив, равномерно)
void gcoordinate_distribution(double * x_g)
{
    for (int i = 0; i < amount; ++i)
    {
        x_g[i] = left + i * (right - left) / (amount - 1);
    }
}

//начальное положение частиц газа (в том числе мнимых)
double fill_image_x(double * image_x)
{
    double new_left = left - (right - left);
    double new_right = right + (right + left);

    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        image_x[j] = new_left + j * (new_right - new_left) / (3*(amount - 1));
    }
}

double spline_kernel(double x_a, double x_b)
{
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

double spline_gradient(double x_a, double x_b)
{
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
    //return (x_a - x_b) * 2. / 3. / h / h * result;
}

//масса газа, находящаяся при постоянной плотности из предположения, что масса всех частиц одинакова
double found_flat_gmass(double * image_x_g, double density)
{
    double mass = 0;
    for (int j = 0; j < 3 * amount - 2; ++j)
    {
        mass += spline_kernel(image_x_g[amount], image_x_g[j]);
    }
    return density / mass;
}

//заполнение массива, содержащего массы частиц газа
void fill_gmass(double * gmass, double * x_g, double * image_x_g, double average_grho)
{
    for(int i = 0; i < amount; ++i)
    {
        //gmass[i] = 1/ (double)amount;
        gmass[i] = gdensity_distribution(x_g[i]) / average_grho * found_flat_gmass(image_x_g, average_grho);
    }
}

void fill_image(double * image, double * real)
{
    for (int i = 0; i < amount; ++i)
    {
        image[i] = real[i];
        image[amount - 1 + i] = real[i];
        image[2 * amount - 2 + i] = real[i];
    }
}

//заполнение массива начального распределения плотности газа
void fill_initial_grho(double * grho, double  * image_gmass, double * x_g, double * image_x_g)
{
    for(int i = 0; i < amount; ++i)
    {
        for (int j = 0; j < 3*amount - 2; ++j)
        {
            grho[i] += image_gmass[j] * spline_kernel(x_g[i], image_x_g[j]);
        }
    }
}

double fill_initial_velocity(double * gvelocity, double * x_g)
{
    for (int i = 0; i < amount; ++i)
    {
        gvelocity[i] = gvelocity_distribution(x_g[i]);
    }
}

double found_next_grho(double * prev_gvelocity, double * gmass, double * x_g, double prev_grho, int i)
{
    double rho = 0;
    for(int j = 0; j < amount; ++j)
    {
        rho += gmass[j] * spline_kernel(x_g[i], x_g[j]);
        //rho += gmass[j] * (prev_gvelocity[i] - prev_gvelocity[j]) * spline_gradient(x_g[i], x_g[j]);
    }
    return rho;
            //tau * rho + prev_grho;
}

double found_next_image_grho(double * image_gmass, double * x_g, double * image_x_g, int i)
{
    double rho = 0;
    for(int j = 0; j < 3*amount - 2; ++j)
    {
        rho += image_gmass[j] * spline_kernel(x_g[i], image_x_g[j]);
    }
    return rho;
}

double found_next_gvelocity(double * prev_grho, double * gmass, double * x_g, double prev_gvelocity, int i)
{
    double velocity = 0;
    for(int j = 0; j < amount; ++j)
    {
        velocity += gmass[j] * (1. / prev_grho[j] + 1. / prev_grho[i]) * spline_gradient(x_g[i], x_g[j]);
    }
    return -tau * pow(c_s, 2) * velocity + prev_gvelocity;
}

double found_next_image_gvelocity(double prev_grho, double * prev_image_grho, double * image_gmass, double * x_g, double * image_x_g,
                                  double prev_gvelocity, int i)
{
    double velocity = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        velocity += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(x_g[i], image_x_g[j]);
    }
    return -tau * pow(c_s, 2) * velocity + prev_gvelocity;
}

double found_next_coordinate(double prev_x_g, double prev_gvelocity)
{
        return prev_x_g + tau * prev_gvelocity;
}

double found_a(double * prev_image_grho, double prev_grho, double * image_gmass, double * image_x_g, double x_g, int i)
{
    double a = 0;
    for(int j = 0; j < 3 * amount - 2; ++j)
    {
        a += image_gmass[j] * (1. / prev_image_grho[j] + 1. / prev_grho) * spline_gradient(x_g, image_x_g[j]);
    }
    return a;

}

int main()
{
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

    gcoordinate_distribution(prev_x_g);
    fill_image_x(prev_image_x_g);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho);
    fill_image(image_gmass, gmass);

    fill_initial_grho(prev_grho, image_gmass, prev_x_g, prev_image_x_g);
    fill_initial_velocity(prev_gvelocity, prev_x_g);

    fill_image(prev_image_grho, prev_grho);
    fill_image(prev_image_gvelocity, prev_gvelocity);

    FILE * fout = fopen("/home/calat/CLionProjects/particles/output.txt", "w");

    for(int i = 0; i < amount; ++i)
    {
        fprintf(fout, "%lf %lf\n", prev_x_g[i], prev_grho[i]);
    }

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
            next_grho[i] = found_next_image_grho(image_gmass, prev_x_g, prev_image_x_g, i);
            next_gvelocity[i] = found_next_image_gvelocity(prev_grho[i], prev_image_grho, image_gmass, prev_x_g, prev_image_x_g, prev_gvelocity[i], i);
            next_x_g[i] = found_next_coordinate(prev_x_g[i], prev_gvelocity[i]);
        }

        fill_image(next_image_grho, next_grho);
        fill_image(next_image_gvelocity, next_gvelocity);

        for (int i = 0; i < 3*amount - 2; ++i)
        {
            next_image_x_g[i] = found_next_coordinate(prev_image_x_g[i], prev_image_gvelocity[i]);
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


    return 0;
}