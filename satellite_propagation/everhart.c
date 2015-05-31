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

const long double R_coeff[7] =
        {
                0.166666666666666667L,
                0.083333333333333333L,
                0.050000000000000000L,
                0.033333333333333333L,
                0.023809523809523810L,
                0.017857142857142857L,
                0.013888888888888889L
        };


const long double V_coeff[7] =
        {
                0.500000000000000000L,
                0.333333333333333333L,
                0.250000000000000000L,
                0.200000000000000000L,
                0.166666666666666667L,
                0.142857142857142857L,
                0.125000000000000000L
        };


void calc_F(long double t, long double pos[3], long double acceleration[3])
{
    long double acc[3];

    acceleration[0] = 0; acceleration[1] = 0; acceleration[2] = 0;

    get_acceleration_by_earth(t, pos[0], pos[1], pos[2], acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

//    get_acceleration_by_moon(t, pos[0], pos[1], pos[2], acc);
//    acceleration[0] += acc[0];
//    acceleration[1] += acc[1];
//    acceleration[2] += acc[2];
//
//    get_acceleration_by_sun(t, pos[0], pos[1], pos[2], acc);
//    acceleration[0] += acc[0];
//    acceleration[1] += acc[1];
//    acceleration[2] += acc[2];

    return;
}


void get_c(long double t[8], long double c[8][8])
//{
//    int i, j;
//
//    for (i = 0; i < N; i++)
//        for (j = 0; j < N; j++)
//            c[i][j] = 0;
//
//    for (i = 0; i < N; i++)
//        c[i][i] = 1;
//
//    for (i = 1; i < N; i++)
//        c[i][0] = -t[i] * c[i-1][0];
//
//    for (i = 1; i < N; i++)
//    {
//        for (j = 1; j < i; j++)
//        {
//            c[i][j] = c[i-1][j-1] - t[i]*c[i-1][j];
//        }
//    }
//}
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            c[i][j] = 0;
        }
        c[i][i] = 1;
    }

    for (i = 1; i < N; i++)
    {
        for (j = 1; j < N; j++)
        {
            c[i][j] = c[i-1][j-1] - t[i]*c[i-1][j];
        }
        c[i][0] = -t[i]*c[i-1][0];
    }

    return;
}


void get_A(long double A[7][3], long double a[7][3], long double c[8][8], int curr_substep)
//{
//    int i, j, k;
//
//    for (i = 0; i < curr_substep; i++)
//    {
//        for (j = 0; j < 3; j++)
//        {
//            A[i][j] = a[i][j];
//        }
//    }
//
//    for (j = 0; j < curr_substep; j++)
//    {
//        for (i = 1; i < 7; i++)
//        {
//            for (k = 0; k < 3; k++)
//            {
//                A[j][k] += c[i][j] * a[i][k];
//            }
//        }
//    }
//}
{
    int i, j, k;

    for (i = 0; i < curr_substep; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A[i][j] = a[i][j];
            for (k = 1; k < 7; k++)
            {
                A[i][j] += c[k][i] * a[k][j];
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
            new_coord[i] += R_coeff[j] * A[j][i] * powl(t, j+3);
            new_vel[i] += V_coeff[j] * A[j][i] * powl(t, j+2);
        }
    }
}


void get_alpha(long double F[3], long double F1[3], long double t[8], long double a[7][3], int n)
{
    int i, j;

    for (i = 0; i < 3; i++)
    {
        a[n][i] = (F[i] - F1[i]) / t[n+1];
    }

    for (i = 1; i < n+1; i++)
    {
        for (j = 0; j < 3; j++)
        {
            a[n][j] -= a[i-1][j];
            a[n][j] /= (t[n+1] - t[i]);
        }
    }
}


void everhart(long double utc_in_mjd,
              long double start_pos[3], long double start_vel[3],
              long double final_pos[3], long double fin_vel[3],
              long double a[7][3], int step)
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

    for (i = 1; i < 9; i++)
    {
        calc_new_pos(t[i], start_pos, start_vel, F1, A, new_pos, new_vel);
        calc_F(utc_in_mjd + t[i]/86400.0, new_pos, F);
        get_alpha(F, F1, t, a, i-1);
        get_A(A, a, c, i);
    }

    calc_new_pos(step, start_pos, start_vel, F1, A, final_pos, fin_vel);

    return;
}