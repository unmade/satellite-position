//
// Created by Леша on 29.05.15.
//

#include <math.h>
#include "forces.h"
#include "date_converters.h"
#include "moon.h"
#include "constants.h"

void get_acceleration_by_moonl(long double utc_in_mjd, long double celes_coord[3],
                               long double acceleration[3])
{
    long double tdb;
    long double moon_coord[3];
    long double r, r3, r0, r03;

    tdb = tt_to_tdbl(mjd_to_ttl(utc_in_mjd));

    get_moon_celestial_positionl(tdb, moon_coord);

    r = sqrtl(powl(moon_coord[0], 2) + powl(moon_coord[1], 2) + powl(moon_coord[2], 2));
    r0 = sqrtl(powl(moon_coord[0]- celes_coord[0], 2) 
               + powl(moon_coord[1] - celes_coord[1], 2)
               + powl(moon_coord[2] - celes_coord[2], 2));
    r3 = powl(r, 3);
    r03 = powl(r0, 3);


    acceleration[0] = FM_M * ((moon_coord[0] - celes_coord[0]) / r03) - FM_M * moon_coord[0] / r3;
    acceleration[1] = FM_M * ((moon_coord[1] - celes_coord[1]) / r03) - FM_M * moon_coord[1] / r3;
    acceleration[2] = FM_M * ((moon_coord[2] - celes_coord[2]) / r03) - FM_M * moon_coord[2] / r3;

    return;
}


void get_acceleration_by_moon(double utc_in_mjd, double celes_coord[3], double acceleration[3])
{
    double tdb;
    double moon_coord[3];
    double r, r3, r0, r03;

    tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));

    get_moon_celestial_position(tdb, moon_coord);

    r = sqrt(pow(moon_coord[0], 2) + pow(moon_coord[1], 2) + pow(moon_coord[2], 2));
    r0 = sqrt(pow(moon_coord[0]- celes_coord[0], 2)
               + pow(moon_coord[1] - celes_coord[1], 2)
               + pow(moon_coord[2] - celes_coord[2], 2));
    r3 = pow(r, 3);
    r03 = pow(r0, 3);


    acceleration[0] = FM_M * ((moon_coord[0] - celes_coord[0]) / r03) - FM_M * moon_coord[0] / r3;
    acceleration[1] = FM_M * ((moon_coord[1] - celes_coord[1]) / r03) - FM_M * moon_coord[1] / r3;
    acceleration[2] = FM_M * ((moon_coord[2] - celes_coord[2]) / r03) - FM_M * moon_coord[2] / r3;

    return;
}