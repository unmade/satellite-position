//
// Created by Леша on 25.05.15.
//

#include <math.h>
#include "coordinates_converters.h"
#include "matrix_operations.h"
#include "date_converters/date_converters.h"
#include "precession.h"
#include "nutation.h"
#include "rotation_matrix.h"

void get_mean_equ_to_fixed_matrixl(long double precession_matrix[3][3], long double m_p[3][3])
{
    transposel(m_p, precession_matrix, 3, 3);
    return;
}

void get_mean_equ_to_fixed_matrix(double precession_matrix[3][3], double m_p[3][3])
{
    transpose(m_p, precession_matrix, 3, 3);
    return;
}


void get_fixed_to_true_equ_matrixl(long double precession_matrix[3][3], long double nutation_matrix[3][3],
                                   long double m_np[3][3])
{
    mult_matricesl(nutation_matrix, precession_matrix, m_np);
    return;
}

void get_fixed_to_true_equ_matrix(double precession_matrix[3][3], double nutation_matrix[3][3],
                                  double m_np[3][3])
{
    mult_matrices(nutation_matrix, precession_matrix, m_np);
    return;
}


void get_true_equ_to_fixed_matrixl(long double m_np[3][3], long double m_pn[3][3])
{
    transposel(m_pn, m_np, 3, 3);
    return;
}

void get_true_equ_to_fixed_matrix(double m_np[3][3], double m_pn[3][3])
{
    transpose(m_pn, m_np, 3, 3);
    return;
}


void get_celes_to_terra_matrixl(long double utc_in_mjd, long double delta_ut, long double m_ct[3][3])
{

    long double precession_matrix[3][3], nutation_matrix[3][3],
                earth_rotation_matrix[3][3], temp[3][3];

    long double tdb = tt_to_tdbl(mjd_to_ttl(utc_in_mjd));
    long double gast = mjd_to_gastl(utc_in_mjd, delta_ut);

    get_precession_matrixl(tdb, precession_matrix);
    get_nutation_matrixl(tdb, nutation_matrix);
    get_earth_rotationl(gast, earth_rotation_matrix);

    mult_matricesl(nutation_matrix, precession_matrix, temp);
    mult_matricesl(earth_rotation_matrix, temp, m_ct);

    return;
}

void get_celes_to_terra_matrix(double utc_in_mjd, double delta_ut, double m_ct[3][3])
{
    double precession_matrix[3][3], nutation_matrix[3][3],
            earth_rotation_matrix[3][3], temp[3][3];

    double tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));
    double gast = mjd_to_gast(utc_in_mjd, delta_ut);

    get_precession_matrix(tdb, precession_matrix);
    get_nutation_matrix(tdb, nutation_matrix);
    get_earth_rotation(gast, earth_rotation_matrix);

    mult_matrices(nutation_matrix, precession_matrix, temp);
    mult_matrices(earth_rotation_matrix, temp, m_ct);

    return;
}


void get_terra_to_celes_matrixl(long double m_ct[3][3], long double m_tc[3][3])
{
    transposel(m_tc, m_ct, 3, 3);
    return;
}

void get_terra_to_celes_matrix(double m_ct[3][3], double m_tc[3][3])
{
    transpose(m_tc, m_ct, 3, 3);
    return;
}


void spherical_to_cartesianl(long double l, long double b, long double r,
                             long double coordinates[3])
{
    coordinates[0] = r * cosl(b) * cosl(l);
    coordinates[1] = r * cosl(b) * sinl(l);
    coordinates[2] = r * sinl(b);

    return;
}


void spherical_to_cartesian(double l, double b, double r,
                             double coordinates[3])
{
    coordinates[0] = r * cos(b) * cos(l);
    coordinates[1] = r * cos(b) * sin(l);
    coordinates[2] = r * sin(b);

    return;
}