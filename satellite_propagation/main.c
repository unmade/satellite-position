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


int main() {
    double a = 0.628318530717959;
    double Rx[3][3];

    rotate_by_x(a, Rx);

    int const n = 3;
    printf("Матрица поворота по ОХ:\n");
    print_matrix(n, n, (double *)Rx);

    double RxT[3][3];
    transpose(RxT, Rx, 3, 3);
    printf("\nТранспонированная матрица поворота по ОХ:\n");
    print_matrix(3, 3, (double *)RxT);


    double jd = utc_to_mjd(2036, 2, 29, 5, 45, 0.0);

    printf("jd=%2.12f\n\n", jd);

    int year, month, day, hour, minute;
    double seconds;
    mjd_to_utc(jd, &year, &month, &day, &hour, &minute, &seconds);

    printf("year=%i, month=%i, day=%i, hour=%i, minute=%i, seconds=%f\n\n",
           year, month, day, hour, minute, seconds);


    printf("∆T = %2.10f\n\n", get_deltaT(2006, 7));


    double tt = utc_to_tt(2003, 11, 15, 15, 35, 0.0);
    printf("TT = %2.12f\n\n", tt);


    double tdb = tt_to_tdb(52958.65004843);
    printf("TDB = %2.12f\n\n", tdb);


    double zeta, theta, z;
    get_precession_parameters(59152, &zeta, &theta, &z);
    printf("ζ(t) = %2.12f; θ(t) = %2.12f; z(t) = %2.12f\n\n",
           zeta*RAD_IN_SEC, theta*RAD_IN_SEC, z*RAD_IN_SEC);


    double precession_matrix[3][3];
    get_precession_matrix(tdb, precession_matrix);
    printf("Precession matrix is:\n");
    print_matrix(3, 3, (double *)precession_matrix);


    double eps_mean = get_eps_mean(59000.5);
    printf("\nε(t) = %2.15f\n\n", eps_mean);


    double fund_args[5];
    printf("Фундаментальные аргументы:\n");
    get_fundumental_args(59000.5, fund_args);
    for (int i=0; i < 5; i++)
    {
        printf("%2.15f   ", fund_args[i]);
    }
    printf("\n\n");


    double delta_psi, delta_eps;
    get_nutation_parameters(55000.5, &delta_psi, &delta_eps);
    printf("∆ψ(t) = %2.15f   ∆ε(t) = %2.15f\n\n", delta_psi, delta_eps);


    double nutation_matrix[3][3];
    printf("Матрица нутации:\n");
    get_nutation_matrix(59152, nutation_matrix);
    print_matrix(3, 3, (double *)nutation_matrix);


    double mjd = utc_to_mjd(2011, 3, 7, 2, 45, 0.0);
    double gmst = utc_to_gmst(mjd, 0);
    printf("\nGMST = %2.15f\n\n", gmst);
//    printf("\nTT = %2.15f\n\n", utc_to_tt(2011,3,7,2,45,0.0));


    double gast = utc_to_gast(2011, 3, 7, 2, 45, 0.0, 0);
    printf("\nGAST = %2.15f\n\n", gast);


    double earth_rotation_matrix[3][3];
    get_earth_rotation_matrix(gast, earth_rotation_matrix);
    printf("Матрица вращения Земли:\n");
    print_matrix(3, 3, (double *)earth_rotation_matrix);



    double tdb2 = tt_to_tdb(utc_to_tt(2003, 9, 23, 7, 13, 0));
    double gast2 = utc_to_gast(2003, 9, 23, 7, 13, 0.0, 0);
    double precession_matrix2[3][3], nutation_matrix2[3][3], earth_rotation_matrix2[3][3],
           m_p[3][3], m_np[3][3], m_pn[3][3], m_ct[3][3], m_tc[3][3];
    get_precession_matrix(tdb2, precession_matrix2);
    get_nutation_matrix(tdb2, nutation_matrix2);
    get_earth_rotation_matrix(gast2, earth_rotation_matrix2);


    printf("\nМатрица перехода от средней подвижной экваториальной к небесной:\n");
    mean_equ_to_fixed(precession_matrix2, m_p);
    print_matrix(3, 3, (double *)m_p);

    printf("\nМатрица перехода от небесной к истинной экватриальной:\n");
    fixed_to_true_equ(precession_matrix2, nutation_matrix2, m_np);
    print_matrix(3, 3, (double *)m_np);

    printf("\nМатрица перехода от истинной экватриальной к небесной:\n");
    true_equ_to_fixed(m_np, m_pn);
    print_matrix(3, 3, (double *)m_pn);

    printf("\nМатрица перехода от небесной к земной:\n");
    tdb2 = tt_to_tdb(utc_to_tt(2011, 3, 7, 2, 45, 0.0));
    gast2 = utc_to_gast(2011, 3, 7, 2, 45, 0.0, 0);
    get_precession_matrix(tdb2, precession_matrix2);
    get_nutation_matrix(tdb2, nutation_matrix2);
    get_earth_rotation_matrix(gast2, earth_rotation_matrix2);
    fixed_to_terra(precession_matrix2, nutation_matrix2, earth_rotation_matrix2, m_ct);
    print_matrix(3, 3, (double *)m_ct);

    printf("\nМатрица перехода от земной к небесной:\n");
    terra_to_fixed(m_ct, m_tc);
    print_matrix(3, 3, (double *)m_tc);


    double utc_in_mjd = 59865.00000000;
    double xc = -43203.606331;
    double yc = 932.853125;
    double zc = 105.030427;

    double Fx, Fy, Fz;

    double sec_in_week = 7*24*60*60;
    double c = 1/86400;
    int i;
    for (i=0; i < sec_in_week; i++) {
        calc(xc, yc, zc, utc_in_mjd, &Fx, &Fy, &Fz);
        utc_in_mjd += c;
//        printf("\nFx = %2.15f\nFy = %2.15f\nFz = %2.15f\n\n", Fx * 1e6, Fy * 1e6, Fz * 1e6);
    }
    printf("\nFx = %2.15f\nFy = %2.15f\nFz = %2.15f\n\n", Fx * 1e6, Fy * 1e6, Fz * 1e6);
//    int x;
//    scanf("%d", &x);
    return 0;
}