//
// Created by Леша on 25.05.15.
//

#include "gravitational_potential.h"
#include "../constants.h"
#include "../coordinates_converters.h"
#include "../matrix_operations.h"


void get_z_and_dzl(long double z, long double r, long double Z[13][13], long double dZ[13][13])
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
        }
        for (n = k + 1; n <= N_MAX; n++)
        {
            dZ[n][k] = (2*n - 1) * Z[n-1][k] + dZ[n-2][k];
        }
    }
}


void get_z_and_dz(double z, double r, double Z[13][13], double dZ[13][13])
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


void get_xy_and_dxyl(long double x, long double y, long double r,
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


void get_xy_and_dxy(double x, double y, double r,
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


    R0_r = R0 / r;
    R = (FM / r) * pow(R0_r, 2);
    dR = 3*FM * pow(R0_r, 2);

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
        dR = (n + 2) * FM * pow(R0_r, n+1);
    }

    u[0] += -FM * (x/r3);
    u[1] += -FM * (y/r3);
    u[2] += -FM * (z/r3);

    get_terra_to_celes_matrix(m_ct, m_tc);
    mult_matrix_by_vector(m_tc, u, acceleration);

    return;
}