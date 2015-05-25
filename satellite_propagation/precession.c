#include "precession.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"
#include "constants.h"


void get_precession_parameters(double tdb, double *zeta, double *theta, double *z)
{
//    double de = (2451545.0-51544.5)/36525.0;
//    double de2 = de*de;
//
    double ts = (tdb - 51544.5) / 36525;
//    double ts2 = ts*ts;
//    double ts3 = ts*ts2;
//
//    double r = (2306.2181+1.39656*de-0.000139*de2)*ts;
//    *zeta = r+(0.30188-0.000344*de)*ts2+0.017998*ts3;
//    *theta = (2004.3109-0.85330*de-0.000217*de2)*ts
//          -(0.42665+0.000217*de)*ts2-0.041833*ts3;
//    *z = r+(1.09468+0.000066*de)*ts2+0.018203*ts3;

    *zeta = SEC_IN_RAD * (2306.2181 + (0.30188 + 0.017998 * ts) * ts) * ts;
    *theta = SEC_IN_RAD * (2004.3109 - (0.42665 + 0.041833 * ts) * ts) * ts;
    *z = SEC_IN_RAD * (2306.2181 + (1.09468 + 0.018203 * ts) * ts) * ts;
    return;
}


void get_precession_matrix(double tdb, double precession_matr[3][3])
{
    double zeta, theta, z;
    get_precession_parameters(tdb, &zeta, &theta, &z);
    double Rz_z[3][3], Ry_theta[3][3], Rz_zeta[3][3], R[3][3];
    rotate_by_z(-z, Rz_z);
    rotate_by_y(theta, Ry_theta);
    rotate_by_z(-zeta, Rz_zeta);

    mult_matrix(Ry_theta, Rz_zeta, R);

    mult_matrix(Rz_z, R, precession_matr);

    return;
}