//
// Created by Леша on 25.05.15.
//

#include <stdio.h>
#include "gravitational_potential.h"
#include "../constants.h"
#include "../coordinates_converters.h"
#include "../date_converters/date_converters.h"
#include "../precession.h"
#include "../nutation.h"
#include "../rotation_matrix.h"
#include "../matrix_operations.h"


void get_r_and_dr(long double r, long double res_R[13], long double res_dR[13])
{
    long double R0_r = R0 / r;
    res_R[0] = res_R[1] = res_dR[0] = res_dR[1] = 0;
    res_R[2] = (FM / r) * powl(R0_r, 2);
    res_dR[2] = 3*FM * powl(R0_r, 2);

    int n;
    for (n = 3; n <= N_MAX; n++)
    {
        res_R[n] = res_R[n -1] * R0_r;
        res_dR[n] = (n + 1) * FM * powl(R0_r, n);
    }
}


void get_z_and_dz(long double z, long double r, long double Z[13][13], long double dZ[13][13])
{
    int n, k;
    for (n = 0; n <= N_MAX; n++)
    {
        for (k = 0; k <= N_MAX; k++)
        {
            Z[n][k] = 0;
            dZ[n][k] = 0;
        }
    }
    Z[0][0] = dZ[0][0] = 0;
    Z[1][0] = z / r;
    dZ[1][0] = 1;
    Z[2][0] = 1.5 * powl(z/r, 2) - 0.5;
    dZ[2][0] = 3* (z / r);


    for (n = 3; n <= N_MAX; n++)
    {
//        Z[n][0] = ((2*n - 1)/n)*(z/r)*Z[n-1][0] - ((n - 1) / n) * Z[n-2][0];
        Z[n][0] = (2*n - 1)*z / (n*r) * Z[n-1][0] - (n-1)*Z[n-2][0] / n;
        dZ[n][0] = n * Z[n-1][0] + (z/r) * dZ[n-1][0];
    }

    for (k = 1; k <= N_MAX; k++)
    {
        for (n = 1; n <= N_MAX; n++)
        {
            Z[n][k] = dZ[n][k-1];
        }
        for (n = 1; n <= k; n++)
        {
            dZ[n][k] = 0;
        }
        for (n = k + 1; n <= N_MAX; n++)
        {
            dZ[n][k] = (2*n - 1) * Z[n-1][k] + dZ[n-2][k];
        }
    }
}


void get_xy_and_dxy(long double x, long double y, long double r,
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


void calc(long double xc, long double yc, long double zc,
          long double utc_in_mjd,
          long double *Fx, long double *Fy, long double *Fz)
{
    long double precession_matrix[3][3], nutation_matrix[3][3],
           earth_rotation_matrix[3][3], m_ct[3][3], coord_matrix[3];

    long double tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));
    long double gast = mjd_to_gast(utc_in_mjd, 0.0);

    get_precession_matrix(tdb, precession_matrix);
    get_nutation_matrix(tdb, nutation_matrix);
    get_earth_rotation_matrix(gast, earth_rotation_matrix);
    fixed_to_terra(precession_matrix, nutation_matrix, earth_rotation_matrix, m_ct);

    long double coord[3];
    coord[0] = xc; coord[1] = yc; coord[2] = zc;
    mult_matrix_by_vector(m_ct, coord, coord_matrix);
    long double x = coord_matrix[0],
                y = coord_matrix[1],
                z = coord_matrix[2];

    long double r, r3, drx,  dry,  drz,
                  dxrx, dxry, dxrz,
                  dyrx, dyry, dyrz,
                  dzrx, dzry, dzrz;

    r = sqrtl(x*x + y*y + z*z);
    r3 = powl(r, 3);

    // ∂(1/r) / ∂x          ∂(1/r) / ∂y             ∂(1/r) / ∂z
    drx = -x/r3;            dry = -y/r3;            drz = -z/r3;

    //∂(x/r) / ∂x           ∂(x/r) / ∂y             ∂(x/r) / ∂z
    dxrx = 1/r - (x*x) / r3;  dxry = -(x*y) / r3;       dxrz = -(x*z) / r3;

    // ∂(y/r) / ∂x          ∂(y/r) / ∂y             ∂(y/r) / ∂z
    dyrx = -(x*y) / r3;       dyry = 1/r - (y*y)/r3;    dyrz = -(y*z) / r3;

    // ∂(z/r) / ∂x          ∂(z/r) / ∂y             ∂(z/r) / ∂z
    dzrx = -(x*z) / r3;       dzry = -(y*z) / r3;       dzrz = 1/r - (z*z) / r3;


    long double R[13], dR[13], Z[13][13], dZ[13][13], X[13], dX_xr[13], dX_yr[13],
            Y[13], dY_xr[13], dY_yr[13];
    get_r_and_dr(r, R, dR);
    get_z_and_dz(z, r, Z, dZ);
    get_xy_and_dxy(x, y, r, X, dX_xr, dX_yr, Y, dY_xr, dY_yr);

    int n, k;
    int harm_index = 0;
    long double Q, dQ_xr, dQ_yr, unk_x = 0, unk_y = 0, unk_z = 0;
    for (n = 2; n <= N_MAX; n++)
    {
        for (k = 0; k <= n; k++)
        {
//            printf("n = %d k = %d ind = %d \n harm_ind = %2.20Lf\n\n", n, k, harm_index, ctes[harm_index]);
            if (k == 0)
            {
                dQ_xr = dQ_yr = 0;
            }
            else
            {
                dQ_xr = ctes[harm_index] * dX_xr[k] + stes[harm_index] * dY_xr[k];
                dQ_yr = ctes[harm_index] * dX_yr[k] + stes[harm_index] * dY_yr[k];
            }

            Q = ctes[harm_index] * X[k] + stes[harm_index] * Y[k];

            unk_x += dR[n] * drx * Z[n][k] * Q   +   R[n] * dZ[n][k] * dzrx * Q
                    + R[n] * Z[n][k] * (dQ_xr * dxrx + dQ_yr * dyrx);
            unk_y += dR[n] * dry * Z[n][k] * Q   +   R[n] * dZ[n][k] * dzry * Q
                     + R[n] * Z[n][k] * (dQ_xr * dxry + dQ_yr * dyry);
            unk_z += dR[n] * drz * Z[n][k] * Q   +   R[n] * dZ[n][k] * dzrz * Q
                     + R[n] * Z[n][k] * (dQ_xr * dxrz + dQ_yr * dyrz);

            ++harm_index;
        }
    }

//    *Fx = -FM * (x/r3);// + unk_x;
//    *Fy = -FM * (y/r3);// + unk_y;
//    *Fz = -FM * (z/r3);// + unk_z;
    *Fx = unk_x;
    *Fy = unk_y;
    *Fz = unk_z;
    return;
}