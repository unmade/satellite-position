//
// Created by Леша on 30.05.15.
//

#include <math.h>
#include "everhart.h"
#include "forces/gravitational_potential.h"
#include "forces/moon_forces.h"
#include "forces/sun_forces.h"


const int N = 8;
const long double h[8] =
        {
                0.000000000000000000L,
                0.056262560526922147L,
                0.180240691736892365L,
                0.352624717113169637L,
                0.547153626330555383L,
                0.734210177215410532L,
                0.885320946839095768L,
                0.977520613561287501L
        };

const long double A_coeff[7] =
        {
                1.0L / 6.0L,
                1.0L / 12.0L,
                1.0L / 20.0L,
                1.0L / 30.0L,
                1.0L / 42.0L,
                1.0L / 56.0L,
                1.0L / 72.0L
        };


const long double V_coeff[7] =
        {
                1.0L / 2.0L,
                1.0L / 3.0L,
                1.0L / 4.0L,
                1.0L / 5.0L,
                1.0L / 6.0L,
                1.0L / 7.0L,
                1.0L / 8.0L
        };


void calc_F(long double t, long double pos[3], long double acceleration[3])
{
    long double acc[3];

    acceleration[0] = 0; acceleration[1] = 0; acceleration[2] = 0;

    calc(pos[0], pos[1], pos[2], t, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_moon(t, pos[0], pos[1], pos[2], acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_sun(t, pos[0], pos[1], pos[2], acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    return;
}


void get_c(long double t[8], long double c[8][8])
{
    int i, j;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            c[i][j] = 0;

    for (i = 0; i < N; i++)
        c[i][i] = 1;

    for (i = 1; i < N; i++)
        c[i][0] = -t[i] * c[i-1][0];

    for (i = 0; i < N; i++)
    {
        for (j = 1; j < i; j++)
        {
            c[i][j] = c[i-1][j-1] - t[i]*c[i-1][j];
        }
    }
}


void get_A(long double A[7][3], long double a[7][3], long double c[8][8], int curr_substep)
{
    int i, j, k;

    for (i = 0; i < curr_substep; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A[i][j] = a[i][j];
        }
    }

    for (j = 0; j < curr_substep; j++)
    {
        for (i = 1; i < 7; i++)
        {
            for (k = 0; k < 3; k++)
            {
                A[j][k] += c[i][j] * a[i][k];
            }
        }
    }
}


void calc_new_pos(long double t,
                  long double pos[3], long double vel[3],
                  long double accel[3], long double A[7][3],
                  long double new_coord[3], long double new_vel[3])
{
    int i, j;

    for (i = 0; i < 3; i++)
    {
        new_coord[i] = pos[i] + vel[i]*t + 0.5 * accel[i] * t*t;
        new_vel[i] = vel[i] + accel[i]*t;
        for (j = 0; j < 7; j++)
        {
            new_coord[i] += A_coeff[j] * A[j][i] * powl(t, j+3);
            new_vel[i] += V_coeff[j] * A[j][i] * powl(t, j+2);
        }
    }
}


void everhart(long double utc_in_mjd, long double start_pos[3], long double start_vel[3], long double final_pos[3], long double fin_vel[3], long double a[7][3], int step)
{
    long double t[8], c[8][8];
    long double A[7][3];
    long double new_pos[3], new_vel[3];
    long double F1[3], F[3];
    int i, j;

    for (i = 0; i < N; i++)
    {
        t[i] = h[i] * step;
    }

    get_c(t, c);


    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A[i][j] = 0;
        }
    }

    calc_F(utc_in_mjd + t[0]/86400.0, start_pos, F1);


    // первый шаг
    calc_new_pos(t[1], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[1]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[0][i] = (F[i] - F1[i]) / t[1];
    }
    get_A(A, a, c, 1);


    // второй шаг
    calc_new_pos(t[2], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[2]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[1][i] = ((F[i] - F1[i]) / t[2] - a[0][i]) / (t[2] - t[1]);
    }
    get_A(A, a, c, 2);

    // третий шаг
    calc_new_pos(t[3], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[3]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[2][i] = (((F[i] - F1[i]) / t[3] - a[0][i]) / (t[3] - t[1]) - a[1][i]) / (t[3] - t[2]);
    }
    get_A(A, a, c, 3);

    // четвертый шаг
    calc_new_pos(t[4], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[4]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[3][i] = ((((F[i] - F1[i]) / t[4] - a[0][i]) / (t[4] - t[1]) - a[1][i]) / (t[4] - t[2]) - a[2][i]) / (t[4] - t[3]);
    }
    get_A(A, a, c, 4);

    // пятый шаг
    calc_new_pos(t[5], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[5]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[4][i] = (((((F[i] - F1[i]) / t[5] - a[0][i]) / (t[5] - t[1]) - a[1][i]) / (t[5] - t[2]) - a[2][i]) / (t[5] - t[3]) - a[3][i]) / (t[5] - t[4]);
    }
    get_A(A, a, c, 5);

    // шестой шаг
    calc_new_pos(t[6], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[6]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[5][i] = ((((((F[i] - F1[i]) / t[6] - a[0][i]) / (t[6] - t[1]) - a[1][i]) / (t[6] - t[2]) - a[2][i]) / (t[6] - t[3]) - a[3][i]) / (t[6] - t[4]) - a[4][i]) / (t[6] - t[5]);
    }
    get_A(A, a, c, 6);

    // седьмой шаг
    calc_new_pos(t[7], start_pos, start_vel, F1, A, new_pos, new_vel);
    calc_F(utc_in_mjd + t[7]/86400.0, new_pos, F);
    for (i = 0; i < 3; i++)
    {
        a[6][i] = (((((((F[i] - F1[i]) / t[7] - a[0][i]) / (t[7] - t[1]) - a[1][i]) / (t[7] - t[2]) - a[2][i]) / (t[7] - t[3]) - a[3][i]) / (t[7] - t[4]) - a[4][i]) / (t[7] - t[5]) - a[5][i]) / (t[7] - t[6]);
    }
    get_A(A, a, c, 7);

    calc_new_pos(step, start_pos, start_vel, F1, A, final_pos, fin_vel);

    return;
}