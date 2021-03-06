//
// Created by Леша on 30.05.15.
//

#include <math.h>
#include "constants.h"
#include "forces.h"


const int N = 8;
static const long double h[8] =
        {
                0.000000000000000000L,
                0.0562625605269221464656522L,
                0.1802406917368923649875799L,
                0.3526247171131696373739078L,
                0.5471536263305553830014486L,
                0.7342101772154105315232106L,
                0.8853209468390957680903598L,
                0.9775206135612875018911745L
        };

static const long double R_coeff[7] =
        {
                0.166666666666666667L,
                0.083333333333333333L,
                0.050000000000000000L,
                0.033333333333333333L,
                0.023809523809523810L,
                0.017857142857142857L,
                0.013888888888888889L
        };


static const long double V_coeff[7] =
        {
                0.500000000000000000L,
                0.333333333333333333L,
                0.250000000000000000L,
                0.200000000000000000L,
                0.166666666666666667L,
                0.142857142857142857L,
                0.125000000000000000L
        };


static void getcl(long double t[8], long double c[8][8])
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

static void getc(double t[8], double c[8][8])
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


static void get_Al(long double A[7][3], long double alpha[7][3], long double c[8][8])
{
    int i, j, k;

    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A[i][j] = alpha[i][j];
            for (k = i+1; k < 7; k++)
            {
                A[i][j] += c[k][i] * alpha[k][j];
            }
        }
    }
}

static void get_A(double A[7][3], double alpha[7][3], double c[8][8])
{
    int i, j, k;

    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A[i][j] = alpha[i][j];
            for (k = i+1; k < 7; k++)
            {
                A[i][j] += c[k][i] * alpha[k][j];
            }
        }
    }
}


static void calc_new_posl(long double t,
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

static void calc_new_pos(double t,
                  double pos[3], double vel[3],
                  double accel[3], double A[7][3],
                  double new_coord[3], double new_vel[3])
{
    int i, j;

    for (i = 0; i < 3; i++)
    {
        new_coord[i] = pos[i] + vel[i]*t + 0.5 * accel[i] * t*t;
        new_vel[i] = vel[i] + accel[i]*t;
        for (j = 0; j < 7; j++)
        {
            new_coord[i] += (double)R_coeff[j] * A[j][i] * powl(t, j+3);
            new_vel[i] += (double)V_coeff[j] * A[j][i] * powl(t, j+2);
        }
    }
}


static void get_alphal(long double F[3], long double F1[3], long double t[8], long double a[7][3], int curr_substep)
{
    int i, j;

    for (j = 0; j < 3; j++)
    {
        a[curr_substep-1][j] = (F[j] - F1[j]) / t[curr_substep];
        for (i = 1; i < curr_substep; i++)
        {
            a[curr_substep - 1][j] -= a[i-1][j];
            a[curr_substep - 1][j] /= (t[curr_substep] - t[i]);
        }
    }
}

static void get_alpha(double F[3], double F1[3], double t[8], double a[7][3], int curr_substep)
{
    int i, j;

    for (j = 0; j < 3; j++)
    {
        a[curr_substep-1][j] = (F[j] - F1[j]) / t[curr_substep];
        for (i = 1; i < curr_substep; i++)
        {
            a[curr_substep - 1][j] -= a[i-1][j];
            a[curr_substep - 1][j] /= (t[curr_substep] - t[i]);
        }
    }
}


void everhartl(long double utc_in_mjd, long double step, long double start_pos[3], long double start_vel[3],
               long double final_pos[3], long double fin_vel[3], long double alpha[7][3])
{
    long double t[8], c[8][8];
    long double A[7][3];
    long double new_pos[3], new_vel[3];
    long double F1[3], F[3];
    int i;
    for (i = 0; i < N; i++)
        t[i] = h[i] * step;

    getcl(t, c);

    get_Al(A, alpha, c);
    get_forcesl(utc_in_mjd + t[0] / SEC_IN_DAY, start_pos, start_vel, F1);

//    for (j = 0; j < 3; j++) {
        for (i = 1; i < N; i++)
        {
            calc_new_posl(t[i], start_pos, start_vel, F1, A, new_pos, new_vel);
            get_forcesl(utc_in_mjd + t[i] / SEC_IN_DAY, new_pos, new_vel, F);
            get_alphal(F, F1, t, alpha, i);
            get_Al(A, alpha, c);
        }
//    }

    calc_new_posl(step, start_pos, start_vel, F1, A, final_pos, fin_vel);

    return;
}

void everhart(double utc_in_mjd, double step, double start_pos[3], double start_vel[3],
              double final_pos[3], double fin_vel[3], double alpha[7][3])
{
    double t[8], c[8][8];
    double A[7][3];
    double new_pos[3], new_vel[3];
    double F1[3], F[3];
    int i, j;

    for (i = 0; i < N; i++)
        t[i] = (double)h[i] * step;

    getc(t, c);

    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            A[i][j] = 0;

    get_forces(utc_in_mjd + t[0] / SEC_IN_DAY, start_pos, start_vel, F1);

    for (i = 1; i < N; i++)
    {
        calc_new_pos(t[i], start_pos, start_vel, F1, A, new_pos, new_vel);
        get_forces(utc_in_mjd + t[i] / SEC_IN_DAY, new_pos, new_vel, F);
        get_alpha(F, F1, t, alpha, i);
        get_A(A, alpha, c);
    }

    calc_new_pos(step, start_pos, start_vel, F1, A, final_pos, fin_vel);

    return;
}