#include <stdio.h>
#include <time.h>

#include "date_converters/date_converters.h"
#include "everhart.h"
#include "tests/rotation_matrix_test.h"
#include "tests/matrix_operations_test.h"
#include "tests/date_converters_test.h"
#include "tests/precession_test.h"
#include "tests/nutation_test.h"
#include "tests/coordinates_converters_test.h"
#include "tests/moon_test.h"
#include "tests/sun_test.h"
#include "tests/forces_test.h"
#include "propagate.h"
#include "tests/propagation_test.h"


int main() {


//    test_rotation_matrices();
//    matrix_operations_test();
//    test_date_converters();
//    test_precession();
//    test_nutation();
//    test_coordinate_converters();
//    test_moon();
//    test_sun();
//    test_forces();
    test_propagation();


    long double start_date = utc_to_mjdl(2015, 5, 14, 6, 0, 0);
    long double end_date = utc_to_mjdl(2015, 5, 15, 6, 0, 0.0);
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
//    long double start_date = utc_to_mjd(2003, 7, 21, 1, 43, 28.2080) - 0.125;
////    double end_date = utc_to_mjd(2003, 7, 24, 1, 19, 58.1360) - 0.125;
//    long double end_date = utc_to_mjd(2003, 7, 28, 00, 48, 37.6900) - 0.125;
////    long double end_date = utc_to_mjd(2003, 7, 24, 01, 19, 58.1360) - 0.125;
////    double end_date = utc_to_mjd(2003, 8, 10, 22, 58, 53.076) - 0.125;
//    long double pos[3], fin_pos[3], vel[3], fin_vel[3];
////    double pos_[3], vel_[3], m_ct[3][3], m_tc[3][3];
//    pos[0] = -7025.194987;
//    pos[1] = -3565.559635;
//    pos[2] = 0.0;
//    vel[0] = 1.453492168;
//    vel[1] = -2.862630493;
//    vel[2] = 6.424834199;

    clock_t start = clock(), diff;

    propagatel(62.0L, start_date, end_date, pos, vel, fin_pos, fin_vel);

    printf("\nx = %2.9f\ny = %2.9f\nz = %2.9f\n\n", (double)fin_pos[0], (double)fin_pos[1], (double)fin_pos[2]);

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