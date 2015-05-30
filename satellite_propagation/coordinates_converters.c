//
// Created by Леша on 25.05.15.
//

#include <math.h>
#include "coordinates_converters.h"
#include "matrix_operations.h"

void mean_equ_to_fixed(long double precession_matrix[3][3], long double m_p[3][3])
{
    transpose(m_p, precession_matrix, 3, 3);
    return;
}


void fixed_to_true_equ(long double precession_matrix[3][3], long double nutation_matrix[3][3],
                       long double m_np[3][3])
{
    mult_matrix(nutation_matrix, precession_matrix, m_np);

    return;
}


void true_equ_to_fixed(long double m_np[3][3], long double m_pn[3][3])
{
    transpose(m_pn, m_np, 3, 3);

    return;
}


void fixed_to_terra(long double precession_matrix[3][3], long double nutation_matrix[3][3],
                    long double earth_rotation_matrix[3][3], long double m_ct[3][3])
{
    long double temp[3][3];
    mult_matrix(nutation_matrix, precession_matrix, temp);
    mult_matrix(earth_rotation_matrix, temp, m_ct);

    return;
}


void terra_to_fixed(long double m_ct[3][3], long double m_tc[3][3])
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