#include "four_eq.h"

double found_next_gvelocity()
{
    return 0;
}

double found_next_dvelocity()
{
    return 0;
}

void whole_system(ParticleParams gas_params, ParticleParams dust_params, ProblemParams problem_params)
{
    int gamount = gas_params.amount;
    int damount = dust_params.amount;

    //Блок массивов для газа.BEGIN
    double prev_x_g[gamount];
    double next_x_g[gamount];

    double gmass[gamount];

    double prev_gvelocity[gamount];
    double next_gvelocity[gamount];
    double prev_grho[gamount];
    double next_grho[gamount];

    double prev_image_x_g[3*gamount - 2];
    double next_image_x_g[3*gamount - 2];

    double image_gmass[3*gamount - 2];

    double prev_image_gvelocity[3*gamount - 2];
    double next_image_gvelocity[3*gamount - 2];
    double prev_image_grho[3*gamount - 2];
    double next_image_grho[3*gamount - 2];

    coordinate_distribution(prev_x_g, gas_params);
    fill_image_x(prev_image_x_g, gas_params);

    double average_grho = gdensity_distribution(0);
    fill_gmass(gmass, prev_x_g, prev_image_x_g, average_grho, gas_params, problem_params);
    fill_image(image_gmass, gmass, gas_params);

    fill_initial_rho(prev_grho, image_gmass, prev_x_g, prev_image_x_g, gas_params, problem_params);
    fill_initial_gvelocity(prev_gvelocity, prev_x_g, gas_params);

    fill_image(prev_image_grho, prev_grho, gas_params);
    fill_image(prev_image_gvelocity, prev_gvelocity, gas_params);
    //Блок массивов для газа.END

    //Блок массивов для пыли.BEGIN
    double prev_x_d[damount];
    double next_x_d[damount];

    double dmass[damount];

    double prev_dvelocity[damount];
    double next_dvelocity[damount];
    double prev_drho[damount];
    double next_drho[damount];

    double prev_image_x_d[3 * damount - 2];
    double next_image_x_d[3 * damount - 2];

    double image_dmass[3 * damount - 2];

    double prev_image_dvelocity[3 * damount - 2];
    double next_image_dvelocity[3 * damount - 2];
    double prev_image_drho[3 * damount - 2];
    double next_image_drho[3 * damount - 2];

    coordinate_distribution(prev_x_d, dust_params);
    fill_image_x(prev_image_x_d, dust_params);

    double average_drho = ddensity_distribution(0);
    fill_dmass(dmass, prev_x_d, prev_image_x_d, average_drho, dust_params, problem_params);
    fill_image(image_dmass, dmass, dust_params);

    fill_initial_rho(prev_drho, image_dmass, prev_x_d, prev_image_x_d, dust_params, problem_params);
    fill_initial_dvelocity(prev_dvelocity, prev_x_d, dust_params);

    fill_image(prev_image_drho, prev_drho, dust_params);
    fill_image(prev_image_dvelocity, prev_dvelocity, dust_params);
    //Блок массивов для пыли.END
}