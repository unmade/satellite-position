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

void get_mean_equ_to_fixed_matrix(long double precession_matrix[3][3], long double m_p[3][3])
{
    transpose(m_p, precession_matrix, 3, 3);
    return;
}


void get_fixed_to_true_equ_matrix(long double precession_matrix[3][3], long double nutation_matrix[3][3],
                                  long double m_np[3][3])
{
    mult_matrix(nutation_matrix, precession_matrix, m_np);

    return;
}


void get_true_equ_to_fixed_matrix(long double m_np[3][3], long double m_pn[3][3])
{
    transpose(m_pn, m_np, 3, 3);

    return;
}


void get_fixed_to_terra_matrix(long double utc_in_mjd, long double m_ct[3][3])
{

    long double precession_matrix[3][3], nutation_matrix[3][3],
                earth_rotation_matrix[3][3], temp[3][3];

    long double tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));
    long double gast = mjd_to_gast(utc_in_mjd, 0.0);

    get_precession_matrix(tdb, precession_matrix);
    get_nutation_matrix(tdb, nutation_matrix);
    get_earth_rotation_matrix(gast, earth_rotation_matrix);

    mult_matrix(nutation_matrix, precession_matrix, temp);
    mult_matrix(earth_rotation_matrix, temp, m_ct);

    return;
}


void get_terra_to_fixed_matrix(long double m_ct[3][3], long double m_tc[3][3])
{
    transpose(m_tc, m_ct, 3, 3);

    return;
}


void spherical_to_cartesian(long double l, long double b, long double r,
                            long double coordinates[3])
{
    coordinates[0] = r * cosl(b) * cosl(l);
    coordinates[1] = r * cosl(b) * sinl(l);
    coordinates[2] = r * sinl(b);

    return;
}