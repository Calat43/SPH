#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

const double pi = 3.14159265358;

const double T = 0.5;
const double left = 0;
const double right = 1;
const int amount = 400;
const double h = 0.08;
const double tau = 0.001;
const double c_s = 1;

//начальное распределение плотности газа
double gdensity_distribution(double x)
{
    return 1;
            //sin(2*pi*x)/100. + 1;
}

//начальное распределение скорости газа
double gvelocity_distribution(double x)
{
    return 0;
    //sin(2.*pi*x)/100.;
}

//начальное положение частиц газа (заполняется сразу как массив, равномерно)
void gcoordinate_distribution(double * x_g)
{
    for (int i = 0; i < amount; ++i)
    {
        x_g[i] = left + i * (right - left) / (amount - 1);
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

double spline_gradient1(double x_a, double x_b)
{
    double r = fabs(x_a - x_b);
    double q = r / h;
    double result = 0;

    if (q >= 0 && q <= 1)
    {
        if (x_a >= x_b)
        {
            result = -3. / h / h * (x_a - x_b) + 9. / 4. / pow(h, 3) * pow((x_a - x_b), 2);
        }
        if (x_a < x_b)
        {
            result = -3. / h / h * (x_a - x_b) - 9. / 4. / pow(h, 3) * pow((x_a - x_b), 2);
        }
    }
    if (q > 1 && q <= 2)
    {
        if(x_a >= x_b)
        {
            result = -3.* x_a/h + 3./h/h * (x_a - x_b) - 3./4./pow(h,3) * pow((x_a - x_b), 2);
        }
        if(x_a < x_b)
        {
            result = 3. * x_a / h + 3./h/h * (x_a - x_b) + 3./4./pow(h, 3) * pow((x_a - x_b), h);
        }
    }
    return 2./3./h * result;
}

//масса газа, находящаяся при постоянной плотности из предположения, что масса всех частиц одинакова
double found_flat_gmass(double * x_g, double density)
{
    double mass = 0;
    for (int j = 0; j < amount; ++j)
    {
        mass += spline_kernel(x_g[amount], x_g[j]);
    }
    return density / mass;
}

//заполнение массива, содержащего массы частиц газа
void fill_gmass(double * gmass, double * x_g, double average_grho)
{
    for(int i = 0; i < amount; ++i)
    {
        gmass[i] = 1/ (double)amount;
        //gmass[i] = gdensity_distribution(x_g[i]) / average_grho * found_flat_gmass(x_g, average_grho);
    }
}

//заполнение массива начального распределения плотности газа
void fill_initial_grho(double * grho, double  * gmass, double * x_g)
{
    for(int i = 0; i < amount; ++i)
    {
        for (int j = 0; j < amount; ++j)
        {
            grho[i] += gmass[j] * spline_kernel(x_g[i], x_g[j]);
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

double found_next_gvelocity(double * prev_grho, double * gmass, double * x_g, double prev_gvelocity, int i)
{
    double velocity = 0;
    for(int j = 0; j < amount; ++j)
    {
        velocity += gmass[j] * (1. / prev_grho[j] + 1. / prev_grho[i]) * spline_gradient(x_g[i], x_g[j]);
    }
    return -tau * pow(c_s, 2) * velocity + prev_gvelocity;
}

double found_next_coordinate(double * prev_x_g, double * prev_gvelocity, int i)
{
        return prev_x_g[i] + tau * prev_gvelocity[i];
}

double found_a(double * prev_grho, double * gmass, double * x_g, int i)
{
    double a = 0;
    for(int j = 0; j < amount; ++j)
    {
        a += gmass[j] * (1. / prev_grho[j] + 1. / prev_grho[i]) * spline_gradient(x_g[i], x_g[j]);
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

    gcoordinate_distribution(prev_x_g);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, average_grho);

    fill_initial_grho(prev_grho, gmass, prev_x_g);
    fill_initial_velocity(prev_gvelocity, prev_x_g);

    double a[amount];

    for (int i = 0; i < amount; ++i)
    {
        a[i] = found_a(prev_grho, gmass, prev_x_g, i);
    }


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
            next_grho[i] = found_next_grho(prev_gvelocity, gmass, prev_x_g, prev_grho[i], i);
            next_gvelocity[i] = found_next_gvelocity(prev_grho, gmass, prev_x_g, prev_gvelocity[i], i);
            next_x_g[i] = found_next_coordinate(prev_x_g, prev_gvelocity, i);
        }

        for(int i = 0; i < amount; ++i)
        {
            prev_grho[i] = next_grho[i];
            prev_gvelocity[i] = next_gvelocity[i];
            prev_x_g[i] = next_x_g[i];

        }
    }


    return 0;
}