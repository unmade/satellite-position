#include "precession_matrix.h"
#include "constants.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"

void calc_precession_parameters(double tdb, double *zeta, double *theta, double *z)
{
    double ts = (tdb - 51544.5) / 36525;

    *zeta = SEC_IN_RAD * (2306.2181 + (0.30188 + 0.017998 * ts) * ts) * ts;
    *theta = SEC_IN_RAD * (2004.3109 - (0.42665 + 0.041833 * ts) * ts) * ts;
    *z = SEC_IN_RAD * (2306.2181 + (1.09468 + 0.018203 * ts) * ts) * ts;
    return;
}


void calc_precession_matrix(double zeta, double theta, double z,
                            double precession_matr[3][3])
{
    double Rz_z[3][3], Ry_theta[3][3], Rz_zeta[3][3], R[3][3];
    rotate_by_z(-z, Rz_z);
    rotate_by_y(theta, Ry_theta);
    rotate_by_z(-zeta, Rz_zeta);

    mult_matrices(Ry_theta, Rz_zeta, R);

    mult_matrices(Rz_z, R, precession_matr);

    return;
}