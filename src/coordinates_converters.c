//
// Created by Леша on 25.05.15.
//

#include <math.h>
#include "coordinates_converters.h"
#include "matrix_operations.h"
#include "date_converters.h"
#include "precession.h"
#include "nutation.h"
#include "rotation_matrix.h"
#include "constants.h"

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


void geodesic_to_terrestriall(long double phi, long double lambda, long double h,
                              long double coordinates[3])
{
    long double g, e2;

    e2 = 2*A0 - A0*A0;
    g = R0 / sqrtl(1 - e2*powl(sinl(phi), 2));
    coordinates[0] = (g + h) * cosl(phi) * cosl(lambda);
    coordinates[1] = (g + h) * cosl(phi) * sinl(lambda);
    coordinates[2] = (g * (1 - e2) + h) * sinl(phi);
    
    return;
}


void geodesic_to_terrestrial(double phi, double lambda, double h,
                             double coordinates[3])
{
    double g, e2;

    e2 = 2*(double)A0 - (double)A0*(double)A0;
    g = (double)R0 / sqrt(1 - e2*pow(sin(phi), 2));
    coordinates[0] = (g + h) * cos(phi) * cos(lambda);
    coordinates[1] = (g + h) * cos(phi) * sin(lambda);
    coordinates[2] = (g * (1 - e2) + h) * sin(phi);

    return;
}


void celestial_to_geodesicl(long double utc_in_mjd, long double coordinates[3], long double geodesic[3])
{
    long double m_ct[3][3], terra_coord[3];
    long double lambdas, phis, phic;
    long double p, q, g, hg, e2;
    int i;

    e2 = 2*A0 - A0*A0;

    get_celes_to_terra_matrixl(utc_in_mjd, 0.0, m_ct);
    mult_matrix_by_vectorl(m_ct, coordinates, terra_coord);

    p = sqrtl(terra_coord[0]*terra_coord[0] + terra_coord[1]*terra_coord[1]);
    lambdas = terra_coord[1] / p;

    hg = 0;
    for (i = 0; i < 3; i++)
    {
        q = sqrtl(powl(terra_coord[2]*(1 + hg), 2) + powl(p * (1 - e2 + hg), 2));
        phis = terra_coord[2]*(1+hg)/q;
        phic = p*(1 - e2 + hg) / q;
        g = R0 / sqrtl(1 - e2 * phis*phis);
        hg = p / (g*phic) - 1;
    }

    geodesic[0] = asinl(phis);
    geodesic[1] = asinl(lambdas);
    geodesic[2] = hg*g;

    return;
}


void celestial_to_geodesic(double utc_in_mjd, double coordinates[3], double geodesic[3])
{
    double m_ct[3][3], terra_coord[3];
    double phis, phic, lambdas;
    double p, q, g, hg, e2;
    int i;

    e2 = 2*(double)A0 - (double)A0*(double)A0;

    get_celes_to_terra_matrix(utc_in_mjd, 0.0, m_ct);
    mult_matrix_by_vector(m_ct, coordinates, terra_coord);

    p = sqrt(terra_coord[0]*terra_coord[0] + terra_coord[1]*terra_coord[1]);
    lambdas = terra_coord[1] / p;

    hg = 0;
    for (i = 0; i < 3; i++)
    {
        q = sqrt(pow(terra_coord[2]*(1 + hg), 2) + pow(p * (1 - e2 + hg), 2));
        phis = terra_coord[2]*(1+hg)/q;
        phic = p*(1 - e2 + hg) / q;
        g = (double)R0 / sqrt(1 - e2 * phis*phis);
        hg = p / (g*phic) - 1;
    }

    geodesic[0] = asin(phis);
    geodesic[1] = asin(lambdas);
    geodesic[2] = hg*g;

    return;
}