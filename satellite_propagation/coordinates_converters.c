//
// Created by Леша on 25.05.15.
//

#include "coordinates_converters.h"
#include "matrix_operations.h"

void mean_equ_to_fixed(double precession_matrix[3][3], double m_p[3][3])
{
    transpose(m_p, precession_matrix, 3, 3);
    return;
}


void fixed_to_true_equ(double precession_matrix[3][3], double nutation_matrix[3][3],
                       double m_np[3][3])
{
    mult_matrix(nutation_matrix, precession_matrix, m_np);
}


void true_equ_to_fixed(double m_np[3][3], double m_pn[3][3])
{
    transpose(m_pn, m_np, 3, 3);
}


void fixed_to_terra(double precession_matrix[3][3], double nutation_matrix[3][3],
                    double earth_rotation_matrix[3][3], double m_ct[3][3])
{
    double temp[3][3];
    mult_matrix(nutation_matrix, precession_matrix, temp);
    mult_matrix(earth_rotation_matrix, temp, m_ct);
}


void terra_to_fixed(double m_ct[3][3], double m_tc[3][3])
{
    transpose(m_tc, m_ct, 3, 3);
}