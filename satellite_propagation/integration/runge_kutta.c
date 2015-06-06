//
// Created by Леша on 05.06.15.
//

#include "integration.h"
#include "../constants.h"
#include "../forces/forces.h"


static const double COEFF[4] = {
    0.0, 0.5, 0.5, 1.0
};


void rungekuttal(long double mjd, long double step,
                 long double r[3], long double v[3],
                 long double fin_r[3], long double fin_v[3])
{
    long double kr[4][3];
    long double kv[4][3];
    long double temp_r[3];
    long double temp_v[3];

    int i, j;

    kr[0][0] = v[0];
    kr[0][1] = v[1];
    kr[0][2] = v[2];

    get_forcesl(mjd, r, v, kv[0]);

    for (j = 1; j < 4; j++)
    {
        for (i = 0; i < 3; i++)
        {
            temp_r[i] = r[i] + COEFF[j] * step * kr[j-1][i];
            temp_v[i] = kr[j][i] = v[i] + COEFF[j] * step * kv[j-1][i];
        }
        get_forcesl(mjd + (step * COEFF[j]) / SEC_IN_DAY, temp_r, temp_v, kv[j]);

    }

    for (i = 0; i < 3; i++)
    {
        fin_r[i] = r[i] + (step / 6) * (kr[0][i] + 2 * kr[1][i] + 2 * kr[2][i] + kr[3][i]);
        fin_v[i] = v[i] + (step / 6) * (kv[0][i] + 2 * kv[1][i] + 2 * kv[2][i] + kv[3][i]);
    }

    return;
}


void rungekutta(double mjd, double step,
                 double r[3], double v[3],
                 double fin_r[3], double fin_v[3])
{
    double kr[4][3];
    double kv[4][3];
    double temp_r[3];
    double temp_v[3];

    int i, j;

    kr[0][0] = v[0];
    kr[0][1] = v[1];
    kr[0][2] = v[2];

    get_forces(mjd, r, v, kv[0]);

    for (j = 1; j < 4; j++)
    {
        for (i = 0; i < 3; i++)
        {
            temp_r[i] = r[i] + COEFF[j] * step * kr[j-1][i];
            temp_v[i] = kr[j][i] = v[i] + COEFF[j] * step * kv[j-1][i];
        }
        get_forces(mjd + (step * COEFF[j]) / SEC_IN_DAY, temp_r, temp_v, kv[j]);

    }

    for (i = 0; i < 3; i++)
    {
        fin_r[i] = r[i] + (step / 6) * (kr[0][i] + 2 * kr[1][i] + 2 * kr[2][i] + kr[3][i]);
        fin_v[i] = v[i] + (step / 6) * (kv[0][i] + 2 * kv[1][i] + 2 * kv[2][i] + kv[3][i]);
    }

    return;
}