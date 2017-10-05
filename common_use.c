#include "common_use.h"

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