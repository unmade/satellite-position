#include "precession.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"
#include "constants.h"


void get_precession_parametersl(long double tdb, long double *zeta, long double *theta, long double *z)
{
    long double ts = (tdb - MJD2000) / JULIAN_C;

    *zeta = SEC_IN_RAD * (2306.2181L + (0.30188L + 0.017998L * ts) * ts) * ts;
    *theta = SEC_IN_RAD * (2004.3109L - (0.42665L + 0.041833L * ts) * ts) * ts;
    *z = SEC_IN_RAD * (2306.2181L + (1.09468L + 0.018203L * ts) * ts) * ts;

    return;
}


void get_precession_parameters(double tdb, double *zeta, double *theta, double *z)
{
    double ts = (tdb - MJD2000) / JULIAN_C;

    *zeta = SEC_IN_RAD * (2306.2181 + (0.30188 + 0.017998 * ts) * ts) * ts;
    *theta = SEC_IN_RAD * (2004.3109 - (0.42665 + 0.041833 * ts) * ts) * ts;
    *z = SEC_IN_RAD * (2306.2181 + (1.09468 + 0.018203 * ts) * ts) * ts;

    return;
}


void get_precession_matrixl(long double tdb, long double precession_matr[3][3])
{
    long double zeta, theta, z;
    long double Rz_z[3][3], Ry_theta[3][3], Rz_zeta[3][3], R[3][3];

    get_precession_parametersl(tdb, &zeta, &theta, &z);

    rotzl(-z, Rz_z);
    rotyl(theta, Ry_theta);
    rotzl(-zeta, Rz_zeta);

    mult_matricesl(Ry_theta, Rz_zeta, R);

    mult_matricesl(Rz_z, R, precession_matr);

    return;
}


void get_precession_matrix(double tdb, double precession_matr[3][3])
{
    double zeta, theta, z;
    double Rz_z[3][3], Ry_theta[3][3], Rz_zeta[3][3], R[3][3];

    get_precession_parameters(tdb, &zeta, &theta, &z);

    rotz(-z, Rz_z);
    roty(theta, Ry_theta);
    rotz(-zeta, Rz_zeta);

    mult_matrices(Ry_theta, Rz_zeta, R);

    mult_matrices(Rz_z, R, precession_matr);

    return;
}