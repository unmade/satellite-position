#include <stdio.h>
#include <time.h>

#include "date_converters.h"
#include "propagate.h"

//#include "rotation_matrix_test.h"
//#include "matrix_operations_test.h"
//#include "date_converters_test.h"
//#include "precession_test.h"
//#include "nutation_test.h"
//#include "coordinates_converters_test.h"
//#include "moon_test.h"
//#include "sun_test.h"
//#include "forces_test.h"
//#include "propagation_test.h"

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
//    test_propagation();


    long double start_date = utc_to_mjd(2015, 5, 14, 6, 0, 0);
    long double end_date = utc_to_mjd(2015, 5, 21, 6, 0, 0.0);
    long double pos[3], fin_pos[3], vel[3], fin_vel[3];

    pos[0] = 33508.4071859207L;
    pos[1] = 25585.7790086467L;
    pos[2] = -493.907358861003L;
    vel[0] = -1.86508722606082L;
    vel[1] = 2.44434537660644L;
    vel[2] = 0.0378635308063818L;

    clock_t start = clock(), diff;


    propagatel(0, 62.0L, start_date, end_date, pos, vel, fin_pos, fin_vel);

    printf("\nx = %2.9f\ny = %2.9f\nz = %2.9f\n\n", (double)fin_pos[0], (double)fin_pos[1], (double)fin_pos[2]);
    printf("vx = %2.9f\nvy = %2.9f\nvz = %2.9f\n\n", (double)fin_vel[0], (double)fin_vel[1], (double)fin_vel[2]);

    diff = clock() - start;

    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);

    return 0;
}