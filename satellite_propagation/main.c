#include <stdio.h>
//#include <math.h>

#include "rotation_matrix.h"
#include "date_converters/date_converters.h"
#include "precession.h"
#include "constants.h"
#include "matrix_operations.h"
#include "nutation.h"
#include "coordinates_converters.h"
#include "forces/gravitational_potential.h"
#include "moon.h"


int main() {
    long double a = 0.628318530717959;
    long double Rx[3][3];

    rotate_by_x(a, Rx);

    int const n = 3;
    printf("Матрица поворота по ОХ:\n");
    print_matrix(n, n, (long double *)Rx);

    long double RxT[3][3];
    transpose(RxT, Rx, 3, 3);
    printf("\nТранспонированная матрица поворота по ОХ:\n");
    print_matrix(3, 3, (long double *)RxT);


    long double jd = utc_to_mjd(2036, 2, 29, 5, 45, 0.0);

    printf("jd=%2.12Lf\n\n", jd);

    int year, month, day, hour, minute;
    long double seconds;
    mjd_to_utc(jd, &year, &month, &day, &hour, &minute, &seconds);

    printf("year=%i, month=%i, day=%i, hour=%i, minute=%i, seconds=%Lf\n\n",
           year, month, day, hour, minute, seconds);


    printf("∆T = %2.10Lf\n\n", get_deltaT(2006, 7));


    long double tt = utc_to_tt(2003, 11, 15, 15, 35, 0.0);
    printf("TT = %2.8Lf\n\n", tt);


    long double tdb = tt_to_tdb(tt);
    printf("TDB = %2.8Lf\n\n", tdb);


    long double zeta, theta, z;
    long double tdb_prec = tt_to_tdb(utc_to_tt(2003, 9, 23, 7, 13, 0.0));
    get_precession_parameters(tdb_prec, &zeta, &theta, &z);
    printf("ζ(t) = %2.5Lf; θ(t) = %2.5Lf; z(t) = %2.5Lf\n\n",
           zeta*RAD_IN_SEC, theta*RAD_IN_SEC, z*RAD_IN_SEC);


    long double precession_matrix[3][3];
    get_precession_matrix(tdb_prec, precession_matrix);
    printf("Precession matrix is:\n");
    print_matrix(3, 3, (long double *)precession_matrix);


    long double eps_mean = get_eps_mean(59000.5);
    printf("\nε(t) = %2.15Lf\n\n", eps_mean);
//
//
    long double fund_args[5];
    int i;
    printf("Фундаментальные аргументы:\n");
    get_fundumental_args(59000.5, fund_args);
    for (i = 0; i < 5; i++)
    {
        printf("%2.15Lf   ", fund_args[i]);
    }
    printf("\n\n");


    long double delta_psi, delta_eps;
    get_nutation_parameters(55000.5, &delta_psi, &delta_eps);
    printf("∆ψ(t) = %2.15Lf   ∆ε(t) = %2.15Lf\n\n", delta_psi, delta_eps);


    long double nutation_matrix[3][3];
    printf("Матрица нутации:\n");
    get_nutation_matrix(55000.5, nutation_matrix);
    print_matrix(3, 3, (long double *)nutation_matrix);

//
//    long double mjd = utc_to_mjd(2011, 3, 7, 2, 45, 0.0);
//    long double gmst = utc_to_gmst(mjd, 0);
//    printf("\nGMST = %2.15Lf\n\n", gmst);
////    printf("\nTT = %2.15f\n\n", utc_to_tt(2011,3,7,2,45,0.0));
//
//
//    long double gast = utc_to_gast(2011, 3, 7, 2, 45, 0.0, 0);
//    printf("\nGAST = %2.15Lf\n\n", gast);
//
//
//    long double earth_rotation_matrix[3][3];
//    get_earth_rotation_matrix(gast, earth_rotation_matrix);
//    printf("Матрица вращения Земли:\n");
//    print_matrix(3, 3, (long double *)earth_rotation_matrix);
//
//
//
//    long double tdb2 = tt_to_tdb(utc_to_tt(2003, 9, 23, 7, 13, 0));
//    long double gast2 = utc_to_gast(2003, 9, 23, 7, 13, 0.0, 0);
//    long double precession_matrix2[3][3], nutation_matrix2[3][3], earth_rotation_matrix2[3][3],
//           m_p[3][3], m_np[3][3], m_pn[3][3], m_ct[3][3], m_tc[3][3];
//    get_precession_matrix(tdb2, precession_matrix2);
//    get_nutation_matrix(tdb2, nutation_matrix2);
//    get_earth_rotation_matrix(gast2, earth_rotation_matrix2);
//
//
//    printf("\nМатрица перехода от средней подвижной экваториальной к небесной:\n");
//    mean_equ_to_fixed(precession_matrix2, m_p);
//    print_matrix(3, 3, (long double *)m_p);
//
//    printf("\nМатрица перехода от небесной к истинной экватриальной:\n");
//    fixed_to_true_equ(precession_matrix2, nutation_matrix2, m_np);
//    print_matrix(3, 3, (long double *)m_np);
//
//    printf("\nМатрица перехода от истинной экватриальной к небесной:\n");
//    true_equ_to_fixed(m_np, m_pn);
//    print_matrix(3, 3, (long double *)m_pn);
//
//    printf("\nМатрица перехода от небесной к земной:\n");
//    tdb2 = tt_to_tdb(utc_to_tt(2011, 3, 7, 2, 45, 0.0));
//    gast2 = utc_to_gast(2011, 3, 7, 2, 45, 0.0, 0);
//    get_precession_matrix(tdb2, precession_matrix2);
//    get_nutation_matrix(tdb2, nutation_matrix2);
//    get_earth_rotation_matrix(gast2, earth_rotation_matrix2);
//    fixed_to_terra(precession_matrix2, nutation_matrix2, earth_rotation_matrix2, m_ct);
//    print_matrix(3, 3, (long double *)m_ct);
//
//    printf("\nМатрица перехода от земной к небесной:\n");
//    terra_to_fixed(m_ct, m_tc);
//    print_matrix(3, 3, (long double *)m_tc);

    printf("Положение луны в эклиптических координатах\n");
//    long double mjd = 58651.11458333;
    long double l, b, r;
    get_moon_ecliptic_position(58651.11536093, &l, &b, &r);
    printf("l = %2.9Lf b = %2.9Lf, r = %2.3Lf", l, b ,r);




//
////
//    long double utc_in_mjd = 59865.00000000L;
//    long double xc = -43203.606331L;
//    long double yc = 932.853125L;
//    long double zc = 105.030427L;

    long double utc_in_mjd = 55433.85;
    long double xc = 4929.325940;
    long double yc = -14681.359180;
    long double zc = 28808.336859;

    long double Fx, Fy, Fz;

//    long double sec_in_week = 7*24*60*60;
//    long double c = 1/86400;
//    int i;
//    for (i=0; i < sec_in_week; i++) {
        calc(xc, yc, zc, utc_in_mjd, &Fx, &Fy, &Fz);
//        utc_in_mjd += c;
//        printf("\nFx = %2.15f\nFy = %2.15f\nFz = %2.15f\n\n", Fx * 1e6, Fy * 1e6, Fz * 1e6);
//    }
    printf("\nFx = %2.9Lf\nFy = %2.9Lf\nFz = %2.9Lf\n\n", Fx * 1e6, Fy * 1e6, Fz * 1e6);
//    int x;
//    scanf("%d", &x);

    return 0;
}