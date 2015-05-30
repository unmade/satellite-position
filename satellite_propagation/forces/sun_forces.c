//
// Created by Леша on 29.05.15.
//

#include "sun_forces.h"
#include "../constants.h"
#include "../sun.h"
#include "../date_converters/date_converters.h"

void get_acceleration_by_sun(long double utc_in_mjd, long double x, long double y, long double z,
                             long double acceleration[3])
{
    long double tdb;
    long double sun_coord[3];
    long double r, r3, r0, r03;

    tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));

    get_sun_celestial_position(tdb, sun_coord);

    sun_coord[0] *= AU;
    sun_coord[1] *= AU;
    sun_coord[2] *= AU;

    r = sqrtl(powl(sun_coord[0], 2) + powl(sun_coord[1], 2) + powl(sun_coord[2], 2));
    r3 = powl(r, 3);
    r0 = sqrtl(powl(sun_coord[0]- x, 2) + powl(sun_coord[1] - y, 2) + powl(sun_coord[2] - z, 2));
    r03 = powl(r0, 3);


    acceleration[0] = FM_S * ((sun_coord[0] - x) / r03) - FM_S * sun_coord[0] / r3;
    acceleration[1] = FM_S * ((sun_coord[1] - y) / r03) - FM_S * sun_coord[1] / r3;
    acceleration[2] = FM_S * ((sun_coord[2] - z) / r03) - FM_S * sun_coord[2] / r3;

    return;
}