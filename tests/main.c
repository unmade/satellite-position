#include <stdio.h>
#include <time.h>

#include "date_converters.h"
#include "propagate.h"

#include "rotation_matrix_test.h"
#include "matrix_operations_test.h"
#include "date_converters_test.h"
#include "precession_test.h"
#include "nutation_test.h"
#include "coordinates_converters_test.h"
#include "moon_test.h"
#include "sun_test.h"
#include "forces_test.h"

int main() {


    test_rotation_matrices();
    matrix_operations_test();
    test_date_converters();
    test_precession();
    test_nutation();
    test_coordinate_converters();
    test_moon();
    test_sun();
    test_forces();
//    test_propagation();


    long double start_date = utc_to_mjd(2015, 5, 14, 6, 0, 0);
//    long double end_date = utc_to_mjdl(2015, 5, 15, 6, 0, 0.0);
    long double end_date = utc_to_mjd(2015, 5, 21, 6, 0, 0.0);
    long double pos[3], fin_pos[3], vel[3], fin_vel[3];

    pos[0] = 33508.4071859207L;
    pos[1] = 25585.7790086467L;
    pos[2] = -493.907358861003L;
    vel[0] = -1.86508722606082L;
    vel[1] = 2.44434537660644L;
    vel[2] = 0.0378635308063818L;

//    pos[0] = 33066.769531250;
//    pos[1] = 26154.355468750;
//    pos[2] = -484.732604980;
//    vel[0] = -1.906557322;
//    vel[1] = 2.412147284;
//    vel[2] = 0.038580179;
//
//    printf("%Lf", (end_date - start_date) * 86400.0);
//
//
    start_date = utc_to_mjd(2003, 7, 21, 1, 43, 28.2080) - 0.125;
//    end_date = utc_to_mjdl(2003, 7, 21, 1, 43, 29.6900) - 0.125;
    end_date = utc_to_mjdl(2003, 7, 28, 00, 48, 37.6900) - 0.125;
//    long double pos[3], fin_pos[3], vel[3], fin_vel[3];
    pos[0] = -7025.194987;
    pos[1] = -3565.559635;
    pos[2] = 0.0;
    vel[0] = 1.453492168;
    vel[1] = -2.862630493;
    vel[2] = 6.424834199;

//    /* Локатор на 14.07.2003 02:52:16.0500 */
//    start_date = utc_to_mjdl(2003, 7, 14, 2, 52, 16.0500L) - 0.125L;
//    pos[0] = -5.964499983e3L;
//    pos[1] = -5.147380808e3L;
//    pos[2] = 0.000000000L;
//    vel[0] = 2.096542117L;
//    vel[1] = -2.428689369L;
//    vel[2] = 6.418883355L;

//    /* Локатор на 17.07.2003 20:08:55.9940 */
//    end_date = utc_to_mjdl(2003, 7, 17, 20, 8, 55.9940) - 0.125L;
//      start_date = utc_to_mjdl(2003, 7, 17, 20, 8, 55.6850L) - 0.125L;
//    pos[0] = -6581.406167706;
//    pos[1] = -4330.333455504;
//    pos[2] = 0.002600154;
//
//    vel[0] = 1.763794244;
//    vel[1] = -2.679890983;
//    vel[2] = 6.419242124;
//    pos[0] = -6.581437439e3L;
//    pos[1] = -4.330317503e3L;
//    pos[2] = 0.000000000L;
//    vel[0] = 1.765289392L;
//    vel[1] = -2.682147192L;
//    vel[2] = 6.424554303L;


    /* Локатор на 18.07.2003 02:06:58.3460 */
//    end_date = utc_to_mjdl(2003, 7, 18, 2, 6, 2.2460L) - 0.125L;
//    end_date = utc_to_mjdl(2003, 7, 18, 2, 6, 58.3460L) - 0.125L;
//    pos[0] = -6.618589941e3L;
//    pos[1] = -4.273473600e3L;
//    pos[2] = 0.000000000L;
//    vel[0] = 1.742045309L;
//    vel[1] = -2.696991698L;
//    vel[2] = 6.424627085L;



    clock_t start = clock(), diff;

//    propagatel(62.0L, start_date, middle_date, pos, vel, fin_pos, fin_vel);
//
//    printf("\nx = %2.9f\ny = %2.9f\nz = %2.9f\n\n", (double)fin_pos[0], (double)fin_pos[1], (double)fin_pos[2]);
//    printf("vx = %2.9f\nvy = %2.9f\nvz = %2.9f\n\n", (double)fin_vel[0], (double)fin_vel[1], (double)fin_vel[2]);
//
//    pos[0] = fin_pos[0];
//    pos[1] = fin_pos[1];
//    pos[2] = fin_pos[2];
//    vel[0] = fin_vel[0];
//    vel[1] = fin_vel[1];
//    vel[2] = fin_vel[2];

    propagatel(100.0L, start_date, end_date, pos, vel, fin_pos, fin_vel);

    printf("\nx = %2.9f\ny = %2.9f\nz = %2.9f\n\n", (double)fin_pos[0], (double)fin_pos[1], (double)fin_pos[2]);
    printf("vx = %2.9f\nvy = %2.9f\nvz = %2.9f\n\n", (double)fin_vel[0], (double)fin_vel[1], (double)fin_vel[2]);

    diff = clock() - start;

    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    int year, month, day, hour, minute;
    long double seconds;
    mjd_to_utcl(start_date + 0.125L, &year, &month, &day, &hour, &minute, &seconds);

    printf("\n\nyear=%i, month=%i, day=%i, hour=%i, minute=%i, seconds=%f\n\n",
           year, month, day, hour, minute, (double)seconds);

    return 0;
}