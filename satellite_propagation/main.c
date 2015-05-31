#include <stdio.h>
#include <time.h>
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
#include "sun.h"
#include "forces/moon_forces.h"
#include "forces/sun_forces.h"
#include "everhart.h"


int main() {
//    long double a = 0.628318530717959;
//    long double Rx[3][3];
//
//    rotate_by_x(a, Rx);
//
//    int const n = 3;
//    printf("Матрица поворота по ОХ:\n");
//    print_matrix(n, n, (long double *)Rx);
//
//    long double RxT[3][3];
//    transpose(RxT, Rx, 3, 3);
//    printf("\nТранспонированная матрица поворота по ОХ:\n");
//    print_matrix(3, 3, (long double *)RxT);
//
//
//    long double jd = utc_to_mjd(2036, 2, 29, 5, 45, 0.0);
//
//    printf("jd=%2.12Lf\n\n", jd);
//
//    int year, month, day, hour, minute;
//    long double seconds;
//    mjd_to_utc(jd, &year, &month, &day, &hour, &minute, &seconds);
//
//    printf("year=%i, month=%i, day=%i, hour=%i, minute=%i, seconds=%Lf\n\n",
//           year, month, day, hour, minute, seconds);
//
//
//    printf("∆T = %2.10Lf\n\n", get_deltaT(2006, 7));
//
//
//    long double tt = utc_to_tt(2003, 11, 15, 15, 35, 0.0);
//    printf("TT = %2.8Lf\n\n", tt);
//
//
//    long double tdb = tt_to_tdb(tt);
//    printf("TDB = %2.8Lf\n\n", tdb);
//
//
//    long double zeta, theta, z;
//    long double tdb_prec = tt_to_tdb(utc_to_tt(2003, 9, 23, 7, 13, 0.0));
//    get_precession_parameters(tdb_prec, &zeta, &theta, &z);
//    printf("ζ(t) = %2.5Lf; θ(t) = %2.5Lf; z(t) = %2.5Lf\n\n",
//           zeta*RAD_IN_SEC, theta*RAD_IN_SEC, z*RAD_IN_SEC);
//
//
//    long double precession_matrix[3][3];
//    get_precession_matrix(tdb_prec, precession_matrix);
//    printf("Precession matrix is:\n");
//    print_matrix(3, 3, (long double *)precession_matrix);
//
//
//    long double eps_mean = get_eps_mean(59000.5);
//    printf("\nε(t) = %2.15Lf\n\n", eps_mean);
////
////
//    long double fund_args[5];
//    int i;
//    printf("Фундаментальные аргументы:\n");
//    get_fundumental_args(59000.5, fund_args);
//    for (i = 0; i < 5; i++)
//    {
//        printf("%2.15Lf   ", fund_args[i]);
//    }
//    printf("\n\n");
//
//
//    long double delta_psi, delta_eps;
//    get_nutation_parameters(55000.5, &delta_psi, &delta_eps);
//    printf("∆ψ(t) = %2.15Lf   ∆ε(t) = %2.15Lf\n\n", delta_psi, delta_eps);
//
//
//    long double nutation_matrix[3][3];
//    printf("Матрица нутации:\n");
//    get_nutation_matrix(55000.5, nutation_matrix);
//    print_matrix(3, 3, (long double *)nutation_matrix);

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
//    get_mean_equ_to_fixed_matrix(precession_matrix2, m_p);
//    print_matrix(3, 3, (long double *)m_p);
//
//    printf("\nМатрица перехода от небесной к истинной экватриальной:\n");
//    get_fixed_to_true_equ_matrix(precession_matrix2, nutation_matrix2, m_np);
//    print_matrix(3, 3, (long double *)m_np);
//
//    printf("\nМатрица перехода от истинной экватриальной к небесной:\n");
//    get_true_equ_to_fixed_matrix(m_np, m_pn);
//    print_matrix(3, 3, (long double *)m_pn);
//
//    printf("\nМатрица перехода от небесной к земной:\n");
//    tdb2 = tt_to_tdb(utc_to_tt(2011, 3, 7, 2, 45, 0.0));
//    gast2 = utc_to_gast(2011, 3, 7, 2, 45, 0.0, 0);
//    get_precession_matrix(tdb2, precession_matrix2);
//    get_nutation_matrix(tdb2, nutation_matrix2);
//    get_earth_rotation_matrix(gast2, earth_rotation_matrix2);
//    get_fixed_to_terra_matrix(precession_matrix2, nutation_matrix2, earth_rotation_matrix2, m_ct);
//    print_matrix(3, 3, (long double *)m_ct);
//
//    printf("\nМатрица перехода от земной к небесной:\n");
//    get_terra_to_fixed_matrix(m_ct, m_tc);
//    print_matrix(3, 3, (long double *)m_tc);

//    printf("Положение Луны в эклиптических координатах:\n");
//    long double mjd = 58651.11458333;
//    long double tdb = 58651.11536093;
//    long double l, b, r;
//    get_moon_ecliptic_position(tdb, &l, &b, &r);
//    printf("l = %2.9Lf b = %2.9Lf, r = %2.3Lf\n\n", l, b, r);

//    long double xm[3];
//    printf("Положение Луны в небесной СК:\n");
//    get_moon_celestial_position(tdb, xm);
//    printf("l = %2.3Lf b = %2.3Lf, r = %2.3Lf\n\n", xm[0], xm[1], xm[2]);
//

//    printf("Положение Солнца в эклиптических координатах:\n");
//    get_sun_ecliptic_position(tdb, &l, &b, &r);
//    printf("l = %2.9Lf b = %2.9Lf, r = %2.3Lf\n\n", l, b, r * AU);


//    printf("Положение Солнца в небесной СК:\n");
//    get_sun_celestial_position(tdb, xm);
//    printf("l = %2.3Lf b = %2.3Lf, r = %2.3Lf\n\n", xm[0]*AU, xm[1]*AU, xm[2]*AU);


//
////
//    long double utc_in_mjd = 59865.00000000 ;
//    long double xc = -43203.606331 ;
//    long double yc = 932.853125 ;
//    long double zc = 105.030427 ;
//    long double acceleration[3];
//    long double utc_in_mjd = 55433.85;
//    long double xc = 4929.325940;
//    long double yc = -14681.359180;
//    long double zc = 28808.336859;

//    get_acceleration_by_earth(utc_in_mjd, xc, yc, zc, acceleration);

//    printf("\nFx = %2.9Lf\nFy = %2.9Lf\nFz = %2.9Lf\n\n", acceleration[0] * 1e6,
//           acceleration[1] * 1e6, acceleration[2] * 1e6);

//    get_acceleration_by_moon(utc_in_mjd, xc, yc, zc, acceleration);
//    printf("\nFx = %2.9Lf\nFy = %2.9Lf\nFz = %2.9Lf\n\n", acceleration[0] * 1e6,
//           acceleration[1] * 1e6, acceleration[2] * 1e6);

//    get_acceleration_by_sun(utc_in_mjd, xc, yc, zc, acceleration);
//    printf("\nFx = %2.9Lf\nFy = %2.9Lf\nFz = %2.9Lf\n\n", acceleration[0] * 1e6,
//           acceleration[1] * 1e6, acceleration[2] * 1e6);


    long double start_date = utc_to_mjd(2015, 5, 14, 6, 0, 0);
//    long double end_date = utc_to_mjd(2015, 5, 15, 6, 0, 0);
    long double pos[3], fin_pos[3], vel[3], fin_vel[3];

//    pos[0] = 33508.4071859207L;
//    pos[1] = 25585.7790086467L;
//    pos[2] = -493.907358861003L;
//    vel[0] = -1.86508722606082L;
//    vel[1] = 2.44434537660644L;
//    vel[2] = 0.0378635308063818L;
    start_date = 51865.0L;
    pos[0] = -43203.607000L;
    pos[1] = 932.853000L;
    pos[2] = 105.030000L;
    vel[0] = -0.108342L;
    vel[1] = -2.994434L;
    vel[2] = 0.062441L;
//    pos_[0] = 33508.4071859207L;
//    pos_[1] = 25585.7790086467L;
//    pos_[2] = -493.907358861003L;
//    vel_[0] = -1.86508722606082L;
//    vel_[1] = 2.44434537660644L;
//    vel_[2] = 0.0378635308063818L;


    clock_t start = clock(), diff;
//
    int step = 1;
    int i, j;
//
    long double a[7][3];
    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 3; j++)
        {
            a[i][j] = 0;
        }
    }
    long double add = (long double)step / 86400.0L;
//    while (start_date <= end_date)
    while (i <= 86400)
    {
        everhart(start_date, pos, vel, fin_pos, fin_vel, a, step);

        pos[0] = fin_pos[0];
        pos[1] = fin_pos[1];
        pos[2] = fin_pos[2];
        vel[0] = fin_vel[0];
        vel[1] = fin_vel[1];
        vel[2] = fin_vel[2];
        start_date += add;
        ++i;
    }

    printf("\nx = %2.9Lf\ny = %2.9Lf\nz = %2.9Lf\n\n", pos[0], pos[1], pos[2]);

    diff = clock() - start;
//
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    return 0;
}