//
// Created by Леша on 25.05.15.
//

#include "constants.h"
#include "coordinates_converters.h"
#include "matrix_operations.h"
#include "forces.h"

static long double const N_MAX = 12; ///< Порядок, до которого учитываются гармоники

static long double const ctes[88] =
        {
                -1.082636023e-3L,   // {  2   0 -1082.63602298 } 1
                -2.41399954e-10L,   // {  2   1 } 2
                +1.57453604e-06L,   // {  2   2 } 3

                +2.532435e-06L,    //  {  3   0     2.53243535 } 4
                +2.19279880e-06L,   // {  3   1 } 5
                +3.09016045e-07L,   // {  3   2 } 6
                +1.00558857e-07L,   // {  3   3 } 7

                +1.619331e-06L,    //  {  4   0     1.61933121 } 8
                -5.08725304e-07L,   // {  4   1 } 9
                +7.84122308e-08L,   // {  4   2 } 10
                +5.92157432e-08L,   // {  4   3 } 11
                -3.98239574e-09L,   // {  4   4 } 12

                +2.277161e-07L,    //  {  5   0     0.22771610 } 13
                -5.37165102e-08L,   // {  5   1 } 14
                +1.05590535e-07L,   // {  5   2 } 15
                -1.49261539e-08L,   // {  5   3 } 16
                -2.29791235e-09L,   // {  5   4 } 17
                +4.30476750e-10L,   // {  5   5 } 18

                -5.396485e-07L,    //  {  6   0    -0.53964849 } 19
                -5.98779769e-08L,   // {  6   1 } 20
                +6.01209884e-09L,   // {  6   2 } 21
                +1.18226641e-09L,   // {  6   3 } 22
                -3.26413891e-10L,   // {  6   4 } 23
                -2.15577115e-10L,   // {  6   5 } 24
                +2.21369255e-12L,   // {  6   6 } 25

                +3.513684e-07L,    //  {  7   0     0.35136844 } 26
                +2.05148728e-07L,   // {  7   1 } 27
                +3.28449048e-08L,   // {  7   2 } 28
                +3.52854052e-09L,   // {  7   3 } 29
                -5.85119492e-10L,   // {  7   4 } 30
                +5.81848555e-13L,   // {  7   5 } 31
                -2.49071768e-11L,   // {  7   6 } 32
                +2.55907810e-14L,   // {  7   7 } 33

                +2.025187e-07L,    //  {  8   0     0.20251872 } 34
                +1.60345872e-08L,   // {  8   1 } 35
                +6.57654233e-09L,   // {  8   2 } 36
                -1.94635815e-10L,   // {  8   3 } 37
                -3.18935802e-10L,   // {  8   4 } 38
                -4.61517343e-12L,   // {  8   5 } 39
                -1.83936427e-12L,   // {  8   6 } 40
                +3.42976182e-13L,   // {  8   7 } 41
                -1.58033229e-13L,   // {  8   8 } 42

                +1.193687e-07L,    //  {  9   0     0.11936871 } 43
                +9.24192722e-08L,   // {  9   1 } 44
                +1.56687424e-09L,   // {  9   2 } 45
                -1.21727527e-09L,   // {  9   3 } 46
                -7.01856112e-12L,   // {  9   4 } 47
                -1.66973699e-12L,   // {  9   5 } 48
                +8.29672520e-13L,   // {  9   6 } 49
                -2.25197343e-13L,   // {  9   7 } 50
                +6.14439391e-14L,   // {  9   8 } 51
                -3.67676273e-15L,   // {  9   9 } 52

                +2.480569e-07L,    //  { 10   0     0.24805686 } 53
                +5.17557870e-08L,   // { 10   1 } 54
                -5.56284564e-09L,   // { 10   2 } 55
                -4.19599893e-11L,   // { 10   3 } 56
                -4.96702541e-11L,   // { 10   4 } 57
                -3.07428287e-12L,   // { 10   5 } 58
                -2.59723183e-13L,   // { 10   6 } 59
                +6.90915376e-15L,   // { 10   7 } 60
                +4.63531420e-15L,   // { 10   8 } 61
                +2.33014817e-15L,   // { 10   9 } 62
                +4.17080240e-16L,   // { 10  10 } 63

                -2.405652e-07L,    //  { 11   0    -0.24056521 } 64
                +9.50842760e-09L,   // { 11   1 } 65
                +9.54203028e-10L,   // { 11   2 } 66
                -1.40960772e-10L,   // { 11   3 } 67
                -1.68525661e-11L,   // { 11   4 } 68
                +1.48944116e-12L,   // { 11   5 } 69
                -5.75467116e-15L,   // { 11   6 } 70
                +1.95426202e-15L,   // { 11   7 } 71
                -2.92494873e-16L,   // { 11   8 } 72
                -1.93432044e-16L,   // { 11   9 } 73
                -4.94639649e-17L,   // { 11  10 } 74
                +9.35170551e-18L,   // { 11  11 } 75

                +1.819117e-07L,     //  { 12   0     0.18191170 } 76
                -3.06800094e-08L,   // { 12   1 } 77
                +6.38039765e-10L,   // { 12   2 } 78
                +1.45191793e-10L,   // { 12   3 } 79
                -2.12381469e-11L,   // { 12   4 } 80
                +8.27990199e-13L,   // { 12   5 } 81
                +7.88309139e-15L,   // { 12   6 } 82
                -4.13155736e-15L,   // { 12   7 } 83
                -5.70825414e-16L,   // { 12   8 } 84
                +1.01272849e-16L,   // { 12   9 } 85
                -1.84017258e-18L,   // { 12  10 } 86
                +4.97869995e-19L,   // { 12  11 } 87
                -2.10894892e-20L    // { 12  12 } 88
        };

static long double const stes[88] =
        {
                0.0L,
                +1.54309997e-09L,   // {  2   1 }
                -9.03868073e-07L,   // {  2   2 }

                0.0L,
                +2.68011894e-07L,   // {  3   1 }
                -2.11402398e-07L,   // {  3   2 }
                +1.97201324e-07L,   // {  3   3 }

                0.0L,
                -4.49459935e-07L,   // {  4   1 }
                +1.48155457e-07L,   // {  4   2 }
                -1.20112918e-08L,   // {  4   3 }
                +6.52560581e-09L,   // {  4   4 }

                0.0L,
                -8.06634638e-08L,   // {  5   1 }
                -5.23267240e-08L,   // {  5   2 }
                -7.10087714e-09L,   // {  5   3 }
                +3.87300507e-10L,   // {  5   4 }
                -1.64820395e-09L,   // {  5   5 }

                0.0L,
                +2.11646643e-08L,   // {  6   1 }
                -4.65039481e-08L,   // {  6   2 }
                +1.84313369e-10L,   // {  6   3 }
                -1.78449133e-09L,   // {  6   4 }
                -4.32918170e-10L,   // {  6   5 }
                -5.52771222e-11L,   // {  6   6 }

                0.0L,
                +6.93698935e-08L,   // {  7   1 }
                +9.28231439e-09L,   // {  7   2 }
                -3.06115024e-09L,   // {  7   3 }
                -2.63618222e-10L,   // {  7   4 }
                +6.39725265e-12L,   // {  7   5 }
                +1.05348786e-11L,   // {  7   6 }
                +4.47598342e-13L,   // {  7   7 }

                0.0L,
                +4.01997816e-08L,   // {  8   1 }
                +5.38131641e-09L,   // {  8   2 }
                -8.72351950e-10L,   // {  8   3 }
                +9.11773560e-11L,   // {  8   4 }
                +1.61252083e-11L,   // {  8   5 }
                +8.62774317e-12L,   // {  8   6 }
                +3.81476567e-13L,   // {  8   7 }
                +1.53533814e-13L,   // {  8   8 }

                0.0L,
                +1.42365696e-08L,   // {  9   1 }
                -2.22867920e-09L,   // {  9   2 }
                -5.63392145e-10L,   // {  9   3 }
                +1.71730872e-11L,   // {  9   4 }
                -5.55091854e-12L,   // {  9   5 }
                +2.94031315e-12L,   // {  9   6 }
                -1.84679217e-13L,   // {  9   7 }
                -9.85618446e-16L,   // {  9   8 }
                +7.44103881e-15L,   // {  9   9 }

                0.0L,
                -8.12891488e-08L,   // { 10   1 }
                -3.05712916e-09L,   // { 10   2 }
                -8.98933286e-10L,   // { 10   3 }
                -4.62248271e-11L,   // { 10   4 }
                -3.12226930e-12L,   // { 10   5 }
                -5.51559139e-13L,   // { 10   6 }
                -2.65068061e-15L,   // { 10   7 }
                -1.05284266e-14L,   // { 10   8 }
                -7.01194816e-16L,   // { 10   9 }
                -9.89260955e-17L,   // { 10  10 }

                0.0L,
                -1.64654645e-08L,   // { 11   1 }
                -5.09736032e-09L,   // { 11   2 }
                -6.86352078e-10L,   // { 11   3 }
                -2.67779792e-11L,   // { 11   4 }
                +1.98250517e-12L,   // { 11   5 }
                +1.34623363e-13L,   // { 11   6 }
                -3.72803733e-14L,   // { 11   7 }
                +1.17044830e-15L,   // { 11   8 }
                +2.58524487e-16L,   // { 11   9 }
                -1.73664923e-17L,   // { 11  10 }
                -1.40785570e-17L,   // { 11  11 }

                0.0L,
                -2.37844845e-08L,   // { 12   1 }
                +1.41642228e-09L,   // { 12   2 }
                +9.15457482e-11L,   // { 12   3 }
                +9.17051709e-13L,   // { 12   4 }
                +2.03324862e-13L,   // { 12   5 }
                +9.33540765e-14L,   // { 12   6 }
                +7.89991291e-15L,   // { 12   7 }
                +3.70152251e-16L,   // { 12   8 }
                +6.13664388e-17L,   // { 12   9 }
                +9.24242436e-18L,   // { 12  10 }
                -2.79007835e-19L,   // { 12  11 }
                -9.83829860e-20L    // { 12  12 }
        };


static void get_z_and_dzl(long double z, long double r, long double Z[13][13], long double dZ[13][13])
{
    long double zr;
    int n, k;
    for (n = 0; n <= N_MAX; n++)
    {
        for (k = 0; k <= N_MAX; k++)
        {
            Z[n][k] = 0;
            dZ[n][k] = 0;
        }
    }

    zr = z / r;
    Z[1][0] = zr;
    dZ[1][0] = 1;
    Z[2][0] = 1.5 * powl(zr, 2) - 0.5;
    dZ[2][0] = 3 * (zr);

    for (n = 3; n <= N_MAX; n++)
    {
        Z[n][0] = (2*n - 1)*z / (n*r) * Z[n-1][0] - (n-1)*Z[n-2][0] / n;
        dZ[n][0] = n * Z[n-1][0] + zr * dZ[n-1][0];
    }

    for (k = 1; k <= N_MAX; k++)
    {
        for (n = 1; n <= N_MAX; n++)
        {
            Z[n][k] = dZ[n][k-1];
            if (n >= k+1) {
                dZ[n][k] = (2 * n - 1) * Z[n - 1][k] + dZ[n - 2][k];
            }
        }
    }
}


static void get_z_and_dz(double z, double r, double Z[13][13], double dZ[13][13])
{
    double zr;
    int n, k;
    for (n = 0; n <= N_MAX; n++)
    {
        for (k = 0; k <= N_MAX; k++)
        {
            Z[n][k] = 0;
            dZ[n][k] = 0;
        }
    }

    zr = z / r;
    Z[1][0] = zr;
    dZ[1][0] = 1;
    Z[2][0] = 1.5 * pow(zr, 2) - 0.5;
    dZ[2][0] = 3 * (zr);
    
    for (n = 3; n <= N_MAX; n++)
    {
        Z[n][0] = (2*n - 1)*z / (n*r) * Z[n-1][0] - (n-1)*Z[n-2][0] / n;
        dZ[n][0] = n * Z[n-1][0] + zr * dZ[n-1][0];
    }

    for (k = 1; k <= N_MAX; k++)
    {
        for (n = 1; n <= N_MAX; n++)
        {
            Z[n][k] = dZ[n][k-1];
        }
        for (n = k + 1; n <= N_MAX; n++)
        {
            dZ[n][k] = (2*n - 1) * Z[n-1][k] + dZ[n-2][k];
        }
    }
}


static void get_xy_and_dxyl(long double x, long double y, long double r,
                     long double X[13], long double dX_xr[13], long double dX_yr[13],
                     long double Y[13], long double dY_xr[13], long double dY_yr[13])
{
    X[0] = 1;
    dX_xr[0] = dX_yr[0] = 0;
    Y[0] = 0;
    dY_xr[0] = dY_yr[0] = 0;
    int k;
    long double xr = x / r;
    long double yr = y / r;

    for(k = 1; k <= N_MAX; k++)
    {
        X[k] = X[k-1]*xr - Y[k-1]*yr;
        dX_xr[k] = dX_xr[k-1]*xr + X[k-1] - dY_xr[k-1]*yr;
        dX_yr[k] = dX_yr[k-1]*xr - dY_yr[k-1]*yr - Y[k-1];

        Y[k] = Y[k-1]*xr + X[k-1]*yr;
        dY_xr[k] = dY_xr[k-1]*xr + Y[k-1] + dX_xr[k-1] * yr;
        dY_yr[k] = dY_yr[k-1]*xr + dX_yr[k-1]*yr + X[k-1];
    }
}


static void get_xy_and_dxy(double x, double y, double r,
                    double X[13], double dX_xr[13], double dX_yr[13],
                    double Y[13], double dY_xr[13], double dY_yr[13])
{
    X[0] = 1;
    dX_xr[0] = dX_yr[0] = 0;
    Y[0] = 0;
    dY_xr[0] = dY_yr[0] = 0;
    int k;
    double xr = x / r;
    double yr = y / r;

    for(k = 1; k <= N_MAX; k++)
    {
        X[k] = X[k-1]*xr - Y[k-1]*yr;
        dX_xr[k] = dX_xr[k-1]*xr + X[k-1] - dY_xr[k-1]*yr;
        dX_yr[k] = dX_yr[k-1]*xr - dY_yr[k-1]*yr - Y[k-1];

        Y[k] = Y[k-1]*xr + X[k-1]*yr;
        dY_xr[k] = dY_xr[k-1]*xr + Y[k-1] + dX_xr[k-1] * yr;
        dY_yr[k] = dY_yr[k-1]*xr + dX_yr[k-1]*yr + X[k-1];
    }
}


void get_acceleration_by_earthl(long double utc_in_mjd, long double celes_coord[3], long double acceleration[3])
{
    long double m_ct[3][3], m_tc[3][3],  // матрицы перехода между земной и небесной СК
                terrestrial_coord[3];

    long double x, y, z; // координаты в земной (прямоугольной) СК

    long double r, r1, r3,
                drx,  dry,  drz,
                dxrx, dxry, dxrz,
                dyrx, dyry, dyrz,
                dzrx, dzry, dzrz;

    long double Z[13][13], dZ[13][13],
                X[13], dX_xr[13], dX_yr[13],
                Y[13], dY_xr[13], dY_yr[13],
                Q, dQ_xr, dQ_yr;

    long double R0_r, R, dR;

    int n, k, harm_index;
    long double u[3];

    get_celes_to_terra_matrixl(utc_in_mjd, 0.0L, m_ct);
    mult_matrix_by_vectorl(m_ct, celes_coord, terrestrial_coord);

    x = terrestrial_coord[0],
    y = terrestrial_coord[1],
    z = terrestrial_coord[2];

    r = sqrtl(x*x + y*y + z*z);
    r1 = 1 / r;
    r3 = powl(r, 3);
    
    // ∂(1/r) / ∂x          ∂(1/r) / ∂y             ∂(1/r) / ∂z
    drx = -x/r3;            dry = -y/r3;            drz = -z/r3;
    //∂(x/r) / ∂x           ∂(x/r) / ∂y             ∂(x/r) / ∂z
    dxrx = r1 - x*x/r3;     dxry = -x*y / r3;       dxrz = -x*z / r3;
    // ∂(y/r) / ∂x          ∂(y/r) / ∂y             ∂(y/r) / ∂z
    dyrx = -x*y / r3;       dyry = r1 - y*y/r3;     dyrz = -y*z / r3;
    // ∂(z/r) / ∂x          ∂(z/r) / ∂y             ∂(z/r) / ∂z
    dzrx = -x*z / r3;       dzry = -y*z / r3;       dzrz = r1 - z*z/r3;


    R0_r = R0 / r;
    R = (FM / r) * powl(R0_r, 2);
    dR = 3*FM * powl(R0_r, 2);

    get_z_and_dzl(z, r, Z, dZ);
    get_xy_and_dxyl(x, y, r, X, dX_xr, dX_yr, Y, dY_xr, dY_yr);

    u[0] = u[1] = u[2] = 0;

    harm_index = 0;
    for (n = 2; n <= N_MAX; n++)
    {
        for (k = 0; k <= n; k++)
        {
            dQ_xr = ctes[harm_index] * dX_xr[k] + stes[harm_index] * dY_xr[k];
            dQ_yr = ctes[harm_index] * dX_yr[k] + stes[harm_index] * dY_yr[k];

            Q = ctes[harm_index] * X[k] + stes[harm_index] * Y[k];

            u[0] += dR * drx * Z[n][k] * Q   +   R * dZ[n][k] * dzrx * Q
                    + R * Z[n][k] * (dQ_xr * dxrx + dQ_yr * dyrx);
            u[1] += dR * dry * Z[n][k] * Q   +   R * dZ[n][k] * dzry * Q
                    + R * Z[n][k] * (dQ_xr * dxry + dQ_yr * dyry);
            u[2] += dR * drz * Z[n][k] * Q   +   R * dZ[n][k] * dzrz * Q
                    + R * Z[n][k] * (dQ_xr * dxrz + dQ_yr * dyrz);

            ++harm_index;
        }
        R *=  R0_r;
        dR = (n + 2) * FM * powl(R0_r, n+1);
    }

    u[0] += -FM * (x/r3);
    u[1] += -FM * (y/r3);
    u[2] += -FM * (z/r3);

    get_terra_to_celes_matrixl(m_ct, m_tc);
    mult_matrix_by_vectorl(m_tc, u, acceleration);

    return;
}


void get_acceleration_by_earth(double utc_in_mjd, double celes_coord[3], double acceleration[3])
{
    double m_ct[3][3], m_tc[3][3],  // матрицы перехода между земной и небесной СК
           terrestrial_coord[3];

    double x, y, z; // координаты в земной (прямоугольной) СК

    double r, r1, r3,
            drx,  dry,  drz,
            dxrx, dxry, dxrz,
            dyrx, dyry, dyrz,
            dzrx, dzry, dzrz;

    double Z[13][13], dZ[13][13],
            X[13], dX_xr[13], dX_yr[13],
            Y[13], dY_xr[13], dY_yr[13],
            Q, dQ_xr, dQ_yr;

    double R0_r, R, dR;

    int n, k, harm_index;
    double u[3];

    get_celes_to_terra_matrix(utc_in_mjd, 0.0, m_ct);
    mult_matrix_by_vector(m_ct, celes_coord, terrestrial_coord);

    x = terrestrial_coord[0],
    y = terrestrial_coord[1],
    z = terrestrial_coord[2];

    r = sqrt(x*x + y*y + z*z);
    r1 = 1 / r;
    r3 = pow(r, 3);

    // ∂(1/r) / ∂x          ∂(1/r) / ∂y             ∂(1/r) / ∂z
    drx = -x/r3;            dry = -y/r3;            drz = -z/r3;
    //∂(x/r) / ∂x           ∂(x/r) / ∂y             ∂(x/r) / ∂z
    dxrx = r1 - x*x/r3;     dxry = -x*y / r3;       dxrz = -x*z / r3;
    // ∂(y/r) / ∂x          ∂(y/r) / ∂y             ∂(y/r) / ∂z
    dyrx = -x*y / r3;       dyry = r1 - y*y/r3;     dyrz = -y*z / r3;
    // ∂(z/r) / ∂x          ∂(z/r) / ∂y             ∂(z/r) / ∂z
    dzrx = -x*z / r3;       dzry = -y*z / r3;       dzrz = r1 - z*z/r3;


    R0_r = (double)R0 / r;
    R = ((double)FM / r) * pow(R0_r, 2);
    dR = 3*(double)FM * pow(R0_r, 2);

    get_z_and_dz(z, r, Z, dZ);
    get_xy_and_dxy(x, y, r, X, dX_xr, dX_yr, Y, dY_xr, dY_yr);

    u[0] = u[1] = u[2] = 0;

    harm_index = 0;
    for (n = 2; n <= N_MAX; n++)
    {
        for (k = 0; k <= n; k++)
        {
            dQ_xr = (double)ctes[harm_index] * dX_xr[k] + (double)stes[harm_index] * dY_xr[k];
            dQ_yr = (double)ctes[harm_index] * dX_yr[k] + (double)stes[harm_index] * dY_yr[k];

            Q = (double)ctes[harm_index] * X[k] + (double)stes[harm_index] * Y[k];

            u[0] += dR * drx * Z[n][k] * Q   +   R * dZ[n][k] * dzrx * Q
                    + R * Z[n][k] * (dQ_xr * dxrx + dQ_yr * dyrx);
            u[1] += dR * dry * Z[n][k] * Q   +   R * dZ[n][k] * dzry * Q
                    + R * Z[n][k] * (dQ_xr * dxry + dQ_yr * dyry);
            u[2] += dR * drz * Z[n][k] * Q   +   R * dZ[n][k] * dzrz * Q
                    + R * Z[n][k] * (dQ_xr * dxrz + dQ_yr * dyrz);

            ++harm_index;
        }
        R *=  R0_r;
        dR = (n + 2) * (double)FM * pow(R0_r, n+1);
    }

    u[0] += -FM * (x/r3);
    u[1] += -FM * (y/r3);
    u[2] += -FM * (z/r3);

    get_terra_to_celes_matrix(m_ct, m_tc);
    mult_matrix_by_vector(m_tc, u, acceleration);

    return;
}