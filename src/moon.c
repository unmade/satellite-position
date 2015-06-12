//
// Created by user on 29.05.2015.
//

#include <math.h>
#include "moon.h"
#include "nutation.h"
#include "constants.h"
#include "coordinates_converters.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"
#include "precession.h"


void get_moon_celestial_positionl(long double tdb, long double coordinates[3])
{
    long double l, b, r;
    long double xe[3];
    long double eps;
    long double Rx[3][3];
    long double precession[3][3], m_pt[3][3], R[3][3];

    get_moon_ecliptic_positionl(tdb, &l, &b, &r);
    spherical_to_cartesianl(l, b, r, xe);

    eps = get_eps_meanl(tdb);

    rotxl(-eps, Rx);

    get_precession_matrixl(tdb, precession);
    transposel(m_pt, precession, 3, 3);

    mult_matricesl(m_pt, Rx, R);
    mult_matrix_by_vectorl(R, xe, coordinates);

    return;
}

void get_moon_celestial_position(double tdb, double coordinates[3])
{
    double l, b, r;
    double xe[3];
    double eps;
    double Rx[3][3];
    double precession[3][3], m_pt[3][3], R[3][3];

    get_moon_ecliptic_position(tdb, &l, &b, &r);
    spherical_to_cartesian(l, b, r, xe);

    eps = get_eps_mean(tdb);

    rotx(-eps, Rx);

    get_precession_matrix(tdb, precession);
    transpose(m_pt, precession, 3, 3);

    mult_matrices(m_pt, Rx, R);
    mult_matrix_by_vector(R, xe, coordinates);

    return;
}


void get_moon_ecliptic_positionl(long double tdb, long double *l, long double *b, long double *r)
{
    long double cfa[5]; // скорректированные фундаментальные аргументы
    long double ar, al, ab;

    get_corr_fund_argsl(tdb, cfa);

    ar= //{3422.70}
          +0.260968 * cosl(4*cfa[4])
         +28.233869 * cosl(2*cfa[4])
          +0.043566 * cosl(cfa[1] + 4*cfa[4])
           +3.08589 * cosl(cfa[1] + 2*cfa[4])
        +186.539296 * cosl(cfa[1])
         +34.311569 * cosl(cfa[1] - 2*cfa[4])
           +0.60071 * cosl(cfa[1] - 4*cfa[4])
          -0.300334 * cosl(cfa[2] + 2*cfa[4])
          -0.399822 * cosl(cfa[2])
          +1.916735 * cosl(cfa[2] - 2*cfa[4])
          +0.034671 * cosl(cfa[2] - 4*cfa[4])
          -0.977818 * cosl(cfa[4])
          +0.282799 * cosl(2*cfa[1] + 2*cfa[4])
         +10.165933 * cosl(2*cfa[1])
          -0.304041 * cosl(2*cfa[1] - 2*cfa[4])
          +0.372337 * cosl(2*cfa[1] - 4*cfa[4])
          -0.949147 * cosl(cfa[1] + cfa[2])
          +1.443617 * cosl(cfa[1] + cfa[2] - 2*cfa[4])
          +0.067283 * cosl(cfa[1] + cfa[2] - 4*cfa[4])
          +0.229935 * cosl(cfa[1] - cfa[2] + 2*cfa[4])
          +1.152852 * cosl(cfa[1] - cfa[2])
          -0.225821 * cosl(cfa[1] - cfa[2] - 2*cfa[4])
          -0.008639 * cosl(2*cfa[2])
          +0.091646 * cosl(2*cfa[2] - 2*cfa[4])
          -0.012103 * cosl(2*cfa[3])
          -0.105291 * cosl(2*cfa[3] - 2*cfa[4])
          -0.109456 * cosl(cfa[1] + cfa[4])
          +0.011715 * cosl(cfa[1] - cfa[4])
          -0.038258 * cosl(cfa[1] - 3*cfa[4])
          +0.149444 * cosl(cfa[2] + cfa[4])
          -0.225821 * cosl(cfa[1] - cfa[2] - 2*cfa[4])
          -0.008639 * cosl(2*cfa[2])
          +0.091646 * cosl(2*cfa[2] - 2*cfa[4])
          -0.012103 * cosl(2*cfa[3])
          -0.105291 * cosl(2*cfa[3] - 2*cfa[4])
          -0.109456 * cosl(cfa[1] + cfa[4])
          +0.011715 * cosl(cfa[1] - cfa[4])
          -0.038258 * cosl(cfa[1] - 3*cfa[4])
          -0.118714 * cosl(3*cfa[1] - 2*cfa[4])
          -0.047853 * cosl(cfa[1] - 2*cfa[3] + 2*cfa[4])
          -0.708093 * cosl(cfa[1] - 2*cfa[3])
          -0.048117 * cosl(cfa[1] + cfa[2] + 2*cfa[4])
          +0.621546 * cosl(3*cfa[1])
          -0.103337 * cosl(2*cfa[1] + cfa[2])
          -0.083228 * cosl(cfa[1] + 2*cfa[3] - 2*cfa[4]);

    //  for longitude in arcsec
    al =  +13.90 * sinl(4*cfa[4])
        +2369.92 * sinl(2*cfa[4])
            +1.98 * sinl(cfa[1] + 4*cfa[4])
          +191.96 * sinl(cfa[1] + 2*cfa[4])
         -4586.47 * sinl(cfa[1] - 2*cfa[4])
           -38.43 * sinl(cfa[1] - 4*cfa[4])
           -24.42 * sinl(cfa[2] + 2*cfa[4])
          -668.15 * sinl(cfa[2])
          -165.15 * sinl(cfa[2] - 2*cfa[4])
            -1.88 * sinl(cfa[2] - 4*cfa[4])
          -125.15 * sinl(cfa[4])
           +14.38 * sinl(2*cfa[1] + 2*cfa[4])
          +769.02 * sinl(2*cfa[1])
          -211.66 * sinl(2*cfa[1] - 2*cfa[4])
           -30.77 * sinl(2*cfa[1] - 4*cfa[4])
            -2.92 * sinl(cfa[1] + cfa[2] + 2*cfa[4])
          -109.67 * sinl(cfa[1] + cfa[2])
          -205.96 * sinl(cfa[1] + cfa[2] - 2*cfa[4])
            -4.39 * sinl(cfa[1] + cfa[2] - 4*cfa[4])
           +14.57 * sinl(cfa[1] - cfa[2] + 2*cfa[4])
          +147.69 * sinl(cfa[1] - cfa[2])
           +28.47 * sinl(cfa[1] - cfa[2] - 2*cfa[4])
            -7.49 * sinl(2*cfa[2])
            -8.09 * sinl(2*cfa[2] - 2*cfa[4])
            -5.74 * sinl(2*cfa[3] + 2*cfa[4])
          -411.60 * sinl(2*cfa[3])
           -55.17 * sinl(2*cfa[3] - 2*cfa[4])
            -8.46 * sinl(cfa[1] + cfa[4])
           +18.61 * sinl(cfa[1] - cfa[4])
            +3.21 * sinl(cfa[1] - 3*cfa[4])
           +18.02 * sinl(cfa[2] + cfa[4])
            +0.56 * sinl(cfa[2] - cfa[4])
           +36.12 * sinl(3*cfa[1])
           -13.19 * sinl(3*cfa[1] - 2*cfa[4])
            -7.65 * sinl(2*cfa[1] + cfa[2])
            +9.70 * sinl(2*cfa[1] - cfa[2])
            -2.49 * sinl(2*cfa[1] - cfa[2] - 2*cfa[4])
            -0.99 * sinl(cfa[1] + 2*cfa[3] + 2*cfa[4])
           -45.10 * sinl(cfa[1] + 2*cfa[3])
            -6.38 * sinl(cfa[1] - 2*cfa[3] + 2*cfa[4])
           +39.53 * sinl(cfa[1] - 2*cfa[3])
            +1.75 * sinl(2*cfa[1] - cfa[4])
        +22639.50 * sinl(cfa[1])
            -0.57 * sinl(2*cfa[1] - 6*cfa[4])
            +0.64 * sinl(cfa[1] - cfa[2] - 4*cfa[4])
            +1.06 * sinl(3*cfa[1] + 2*cfa[4])
            -1.19 * sinl(3*cfa[1] - 4*cfa[4])
            -8.63 * sinl(2*cfa[1] + cfa[2] - 2*cfa[4])
            -2.74 * sinl(2*cfa[1] + cfa[2] - 4*cfa[4])
            +1.18 * sinl(2*cfa[1] - cfa[2] + 2*cfa[4])
            -1.17 * sinl(cfa[1] + 2*cfa[2])
            -7.41 * sinl(cfa[1] + 2*cfa[2] - 2*cfa[4])
            +0.76 * sinl(cfa[1] - 2*cfa[2] + 2*cfa[4])
            +2.58 * sinl(cfa[1] - 2*cfa[2])
            +2.53 * sinl(cfa[1] - 2*cfa[2] - 2*cfa[4])
            +9.37 * sinl(cfa[1] - 2*cfa[3] - 2*cfa[4])
            -2.15 * sinl(cfa[2] + 2*cfa[3] - 2*cfa[4])
            -1.44 * sinl(cfa[2] - 2*cfa[3] + 2*cfa[4])
            -0.59 * sinl(2*cfa[1] + cfa[4])
            +1.22 * sinl(2*cfa[1] - 3*cfa[4])
            +1.27 * sinl(cfa[1]+cfa[2] + cfa[4])
            -1.09 * sinl(cfa[1]-cfa[2] - cfa[4])
            +0.58 * sinl(2*cfa[3] - cfa[4])
            +1.94 * sinl(4*cfa[1])
            -0.95 * sinl(4*cfa[1] - 2*cfa[4])
            -0.55 * sinl(3*cfa[1] + cfa[2])
            +0.67 * sinl(3*cfa[1] - cfa[2])
            -4.00 * sinl(2*cfa[1] + 2*cfa[3])
            +0.56 * sinl(2*cfa[1] + 2*cfa[3] - 2*cfa[4])
            -1.30 * sinl(2*cfa[1] - 2*cfa[3])
            +0.54 * sinl(2*cfa[1] - 2*cfa[3] - 2*cfa[4]);

    //  for latitude in arcsec
    ab =     117.26 * sinl(cfa[3] + 2*cfa[4])
          +18461.35 * sinl(cfa[3])
            -623.66 * sinl(cfa[3] - 2*cfa[4])
              -3.67 * sinl(cfa[3] - 4*cfa[4])
             +15.12 * sinl(cfa[1] + cfa[3] + 2*cfa[4])
            -166.58 * sinl(cfa[1] + cfa[3] - 2*cfa[4])
              -6.58 * sinl(cfa[1] + cfa[3] - 4*cfa[4])
              +3.00 * sinl(-cfa[1] + cfa[3] + 4*cfa[4])
            +199.49 * sinl(-cfa[1] + cfa[3] + 2*cfa[4])
            -999.69 * sinl(-cfa[1] + cfa[3])
             -33.36 * sinl(-cfa[1] + cfa[3] - 2*cfa[4])
              -6.48 * sinl(cfa[2] + cfa[3])
             -29.65 * sinl(cfa[2] + cfa[3] - 2*cfa[4])
              +7.98 * sinl(-cfa[2] + cfa[3] + 2*cfa[4])
              +4.86 * sinl(-cfa[2] + cfa[3])
              -5.38 * sinl(cfa[3] + cfa[4])
              +4.81 * sinl(cfa[3] - cfa[4])
             -15.57 * sinl(2*cfa[1] + cfa[3] - 2*cfa[4])
             -31.76 * sinl(-2*cfa[1] + cfa[3])
              -5.33 * sinl(cfa[1] + cfa[2] + cfa[3])
              +8.89 * sinl(-cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
              +6.75 * sinl(cfa[1] - cfa[2] + cfa[3])
              -5.65 * sinl(-cfa[1] + cfa[2] + cfa[3])
              -1.02 * sinl(cfa[1] + 3*cfa[3])
              +1.19 * sinl(cfa[3] + 4*cfa[4])
           +1010.16 * sinl(cfa[1] + cfa[3])
              -1.26 * sinl(cfa[2] + cfa[3] + 2*cfa[4])
             +12.12 * sinl(-cfa[2] + cfa[3] - 2*cfa[4])
              -6.29 * sinl(3*cfa[3])
              -2.18 * sinl(3*cfa[3] - 2*cfa[4])
              +1.51 * sinl(2*cfa[1] + cfa[3] + 2*cfa[4])
             +61.91 * sinl(2*cfa[1] + cfa[3])
              +2.41 * sinl(-2*cfa[1] + cfa[3] + 4*cfa[4])
              -1.62 * sinl(-2*cfa[1] + cfa[3] + 2*cfa[4])
              -2.14 * sinl(-2*cfa[1] + cfa[3] - 2*cfa[4])
              -7.45 * sinl(cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
              +5.08 * sinl(-cfa[1] - cfa[2] + cfa[3])
              +1.13 * sinl(cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
              -1.32 * sinl(-cfa[1] + cfa[2] + cfa[3] + 2*cfa[4])
              -1.77 * sinl(-cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
              -1.09 * sinl(2*cfa[2] + cfa[3] - 2*cfa[4])
              -2.79 * sinl(-cfa[1] + 3*cfa[3])
              +3.98 * sinl(3*cfa[1] + cfa[3])
              -1.51 * sinl(3*cfa[1] + cfa[3] - 2*cfa[4])
              -1.58 * sinl(-3*cfa[1] + cfa[3]);

    *r = R0 / (SEC_IN_RAD*0.999953253 * (3422.7 + ar));
    *l = cfa[0] + SEC_IN_RAD*al;
    *b = SEC_IN_RAD * ab;
}

void get_moon_ecliptic_position(double tdb,double *l, double *b, double *r)
{
    double cfa[5]; // скорректированные фундаментальные аргументы
    double ar, al, ab;

    get_corr_fund_args(tdb, cfa);

    ar= //{3422.70}
              +0.260968 * cos(4*cfa[4])
             +28.233869 * cos(2*cfa[4])
              +0.043566 * cos(cfa[1] + 4*cfa[4])
               +3.08589 * cos(cfa[1] + 2*cfa[4])
            +186.539296 * cos(cfa[1])
             +34.311569 * cos(cfa[1] - 2*cfa[4])
               +0.60071 * cos(cfa[1] - 4*cfa[4])
              -0.300334 * cos(cfa[2] + 2*cfa[4])
              -0.399822 * cos(cfa[2])
              +1.916735 * cos(cfa[2] - 2*cfa[4])
              +0.034671 * cos(cfa[2] - 4*cfa[4])
              -0.977818 * cos(cfa[4])
              +0.282799 * cos(2*cfa[1] + 2*cfa[4])
             +10.165933 * cos(2*cfa[1])
              -0.304041 * cos(2*cfa[1] - 2*cfa[4])
              +0.372337 * cos(2*cfa[1] - 4*cfa[4])
              -0.949147 * cos(cfa[1] + cfa[2])
              +1.443617 * cos(cfa[1] + cfa[2] - 2*cfa[4])
              +0.067283 * cos(cfa[1] + cfa[2] - 4*cfa[4])
              +0.229935 * cos(cfa[1] - cfa[2] + 2*cfa[4])
              +1.152852 * cos(cfa[1] - cfa[2])
              -0.225821 * cos(cfa[1] - cfa[2] - 2*cfa[4])
              -0.008639 * cos(2*cfa[2])
              +0.091646 * cos(2*cfa[2] - 2*cfa[4])
              -0.012103 * cos(2*cfa[3])
              -0.105291 * cos(2*cfa[3] - 2*cfa[4])
              -0.109456 * cos(cfa[1] + cfa[4])
              +0.011715 * cos(cfa[1] - cfa[4])
              -0.038258 * cos(cfa[1] - 3*cfa[4])
              +0.149444 * cos(cfa[2] + cfa[4])
              -0.225821 * cos(cfa[1] - cfa[2] - 2*cfa[4])
              -0.008639 * cos(2*cfa[2])
              +0.091646 * cos(2*cfa[2] - 2*cfa[4])
              -0.012103 * cos(2*cfa[3])
              -0.105291 * cos(2*cfa[3] - 2*cfa[4])
              -0.109456 * cos(cfa[1] + cfa[4])
              +0.011715 * cos(cfa[1] - cfa[4])
              -0.038258 * cos(cfa[1] - 3*cfa[4])
              -0.118714 * cos(3*cfa[1] - 2*cfa[4])
              -0.047853 * cos(cfa[1] - 2*cfa[3] + 2*cfa[4])
              -0.708093 * cos(cfa[1] - 2*cfa[3])
              -0.048117 * cos(cfa[1] + cfa[2] + 2*cfa[4])
              +0.621546 * cos(3*cfa[1])
              -0.103337 * cos(2*cfa[1] + cfa[2])
              -0.083228 * cos(cfa[1] + 2*cfa[3] - 2*cfa[4]);

    //  for longitude in arcsec
    al =  +13.90 * sin(4*cfa[4])
          +2369.92 * sin(2*cfa[4])
          +1.98 * sin(cfa[1] + 4*cfa[4])
          +191.96 * sin(cfa[1] + 2*cfa[4])
          -4586.47 * sin(cfa[1] - 2*cfa[4])
          -38.43 * sin(cfa[1] - 4*cfa[4])
          -24.42 * sin(cfa[2] + 2*cfa[4])
          -668.15 * sin(cfa[2])
          -165.15 * sin(cfa[2] - 2*cfa[4])
          -1.88 * sin(cfa[2] - 4*cfa[4])
          -125.15 * sin(cfa[4])
          +14.38 * sin(2*cfa[1] + 2*cfa[4])
          +769.02 * sin(2*cfa[1])
          -211.66 * sin(2*cfa[1] - 2*cfa[4])
          -30.77 * sin(2*cfa[1] - 4*cfa[4])
          -2.92 * sin(cfa[1] + cfa[2] + 2*cfa[4])
          -109.67 * sin(cfa[1] + cfa[2])
          -205.96 * sin(cfa[1] + cfa[2] - 2*cfa[4])
          -4.39 * sin(cfa[1] + cfa[2] - 4*cfa[4])
          +14.57 * sin(cfa[1] - cfa[2] + 2*cfa[4])
          +147.69 * sin(cfa[1] - cfa[2])
          +28.47 * sin(cfa[1] - cfa[2] - 2*cfa[4])
          -7.49 * sin(2*cfa[2])
          -8.09 * sin(2*cfa[2] - 2*cfa[4])
          -5.74 * sin(2*cfa[3] + 2*cfa[4])
          -411.60 * sin(2*cfa[3])
          -55.17 * sin(2*cfa[3] - 2*cfa[4])
          -8.46 * sin(cfa[1] + cfa[4])
          +18.61 * sin(cfa[1] - cfa[4])
          +3.21 * sin(cfa[1] - 3*cfa[4])
          +18.02 * sin(cfa[2] + cfa[4])
          +0.56 * sin(cfa[2] - cfa[4])
          +36.12 * sin(3*cfa[1])
          -13.19 * sin(3*cfa[1] - 2*cfa[4])
          -7.65 * sin(2*cfa[1] + cfa[2])
          +9.70 * sin(2*cfa[1] - cfa[2])
          -2.49 * sin(2*cfa[1] - cfa[2] - 2*cfa[4])
          -0.99 * sin(cfa[1] + 2*cfa[3] + 2*cfa[4])
          -45.10 * sin(cfa[1] + 2*cfa[3])
          -6.38 * sin(cfa[1] - 2*cfa[3] + 2*cfa[4])
          +39.53 * sin(cfa[1] - 2*cfa[3])
          +1.75 * sin(2*cfa[1] - cfa[4])
          +22639.50 * sin(cfa[1])
          -0.57 * sin(2*cfa[1] - 6*cfa[4])
          +0.64 * sin(cfa[1] - cfa[2] - 4*cfa[4])
          +1.06 * sin(3*cfa[1] + 2*cfa[4])
          -1.19 * sin(3*cfa[1] - 4*cfa[4])
          -8.63 * sin(2*cfa[1] + cfa[2] - 2*cfa[4])
          -2.74 * sin(2*cfa[1] + cfa[2] - 4*cfa[4])
          +1.18 * sin(2*cfa[1] - cfa[2] + 2*cfa[4])
          -1.17 * sin(cfa[1] + 2*cfa[2])
          -7.41 * sin(cfa[1] + 2*cfa[2] - 2*cfa[4])
          +0.76 * sin(cfa[1] - 2*cfa[2] + 2*cfa[4])
          +2.58 * sin(cfa[1] - 2*cfa[2])
          +2.53 * sin(cfa[1] - 2*cfa[2] - 2*cfa[4])
          +9.37 * sin(cfa[1] - 2*cfa[3] - 2*cfa[4])
          -2.15 * sin(cfa[2] + 2*cfa[3] - 2*cfa[4])
          -1.44 * sin(cfa[2] - 2*cfa[3] + 2*cfa[4])
          -0.59 * sin(2*cfa[1] + cfa[4])
          +1.22 * sin(2*cfa[1] - 3*cfa[4])
          +1.27 * sin(cfa[1]+cfa[2] + cfa[4])
          -1.09 * sin(cfa[1]-cfa[2] - cfa[4])
          +0.58 * sin(2*cfa[3] - cfa[4])
          +1.94 * sin(4*cfa[1])
          -0.95 * sin(4*cfa[1] - 2*cfa[4])
          -0.55 * sin(3*cfa[1] + cfa[2])
          +0.67 * sin(3*cfa[1] - cfa[2])
          -4.00 * sin(2*cfa[1] + 2*cfa[3])
          +0.56 * sin(2*cfa[1] + 2*cfa[3] - 2*cfa[4])
          -1.30 * sin(2*cfa[1] - 2*cfa[3])
          +0.54 * sin(2*cfa[1] - 2*cfa[3] - 2*cfa[4]);

    //  for latitude in arcsec
    ab =     117.26 * sin(cfa[3] + 2*cfa[4])
             +18461.35 * sin(cfa[3])
             -623.66 * sin(cfa[3] - 2*cfa[4])
             -3.67 * sin(cfa[3] - 4*cfa[4])
             +15.12 * sin(cfa[1] + cfa[3] + 2*cfa[4])
             -166.58 * sin(cfa[1] + cfa[3] - 2*cfa[4])
             -6.58 * sin(cfa[1] + cfa[3] - 4*cfa[4])
             +3.00 * sin(-cfa[1] + cfa[3] + 4*cfa[4])
             +199.49 * sin(-cfa[1] + cfa[3] + 2*cfa[4])
             -999.69 * sin(-cfa[1] + cfa[3])
             -33.36 * sin(-cfa[1] + cfa[3] - 2*cfa[4])
             -6.48 * sin(cfa[2] + cfa[3])
             -29.65 * sin(cfa[2] + cfa[3] - 2*cfa[4])
             +7.98 * sin(-cfa[2] + cfa[3] + 2*cfa[4])
             +4.86 * sin(-cfa[2] + cfa[3])
             -5.38 * sin(cfa[3] + cfa[4])
             +4.81 * sin(cfa[3] - cfa[4])
             -15.57 * sin(2*cfa[1] + cfa[3] - 2*cfa[4])
             -31.76 * sin(-2*cfa[1] + cfa[3])
             -5.33 * sin(cfa[1] + cfa[2] + cfa[3])
             +8.89 * sin(-cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
             +6.75 * sin(cfa[1] - cfa[2] + cfa[3])
             -5.65 * sin(-cfa[1] + cfa[2] + cfa[3])
             -1.02 * sin(cfa[1] + 3*cfa[3])
             +1.19 * sin(cfa[3] + 4*cfa[4])
             +1010.16 * sin(cfa[1] + cfa[3])
             -1.26 * sin(cfa[2] + cfa[3] + 2*cfa[4])
             +12.12 * sin(-cfa[2] + cfa[3] - 2*cfa[4])
             -6.29 * sin(3*cfa[3])
             -2.18 * sin(3*cfa[3] - 2*cfa[4])
             +1.51 * sin(2*cfa[1] + cfa[3] + 2*cfa[4])
             +61.91 * sin(2*cfa[1] + cfa[3])
             +2.41 * sin(-2*cfa[1] + cfa[3] + 4*cfa[4])
             -1.62 * sin(-2*cfa[1] + cfa[3] + 2*cfa[4])
             -2.14 * sin(-2*cfa[1] + cfa[3] - 2*cfa[4])
             -7.45 * sin(cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
             +5.08 * sin(-cfa[1] - cfa[2] + cfa[3])
             +1.13 * sin(cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
             -1.32 * sin(-cfa[1] + cfa[2] + cfa[3] + 2*cfa[4])
             -1.77 * sin(-cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
             -1.09 * sin(2*cfa[2] + cfa[3] - 2*cfa[4])
             -2.79 * sin(-cfa[1] + 3*cfa[3])
             +3.98 * sin(3*cfa[1] + cfa[3])
             -1.51 * sin(3*cfa[1] + cfa[3] - 2*cfa[4])
             -1.58 * sin(-3*cfa[1] + cfa[3]);

    *r = (double)R0 / (SEC_IN_RAD*0.999953253 * (3422.7 + ar));
    *l = cfa[0] + SEC_IN_RAD*al;
    *b = SEC_IN_RAD * ab;
}
