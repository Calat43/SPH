#include "fscanf.h"

void xy_value(char * file_in, char * file_out)
{
    FILE * dataFile = fopen(file_out, "r");
    double coord[100];
    double vgas[100];
    double vdust[100];
    double rhogas[100];
    double rhodust[100];

    double x[100];
    double y[100];


    char ignore[1024];
    fgets(ignore, sizeof(ignore), dataFile);

    for (int i = 0; i < 100; ++i)
    {
        fscanf(dataFile, "%lf %lf %lf %lf %lf", &(coord[i]), &(vgas[i]), &(vdust[i]), &(rhogas[i]), &(rhodust[i]));
    }

    for (int i = 0; i < 100; ++i)
    {
        x[i] = vgas[i] - vdust[i];
        y[i] = vgas[i] + rhodust[i] / rhogas[i] * vdust[i];
    }

    FILE * fout = fopen(file_in, "w");

    for (int i = 0; i < 100; ++i)
    {
        fprintf(fout, "%0.8lf %0.8lf %0.8lf \n", coord[i], x[i], y[i]);
    }
}

