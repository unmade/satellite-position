//
// Created by Леша on 04.06.15.
//

#include <math.h>
#include "forces.h"
#include "constants.h"
#include "sun.h"
#include "date_converters.h"

void get_acceleration_by_atmospherel(long double coord[3], long double vel[3], long double acceleration[3])
{
    long double rc, p, h, vc;
    int i;
    rc = sqrtl(powl(coord[0], 2) + powl(coord[1], 2) + powl(coord[2], 2));

    h = rc - R0 * (1 - A0*powl(coord[2]/rc, 2));
    p = 2e-13L * expl(-((h-200)/60));

    vc = sqrtl(powl(vel[0], 2) + powl(vel[1], 2) + powl(vel[2], 2));

    for (i = 0; i < 3; i++) {
        acceleration[i] = -Sb * p * vc * vel[i];
    }

    return;
}


void get_acceleration_by_atmosphere(double coord[3], double vel[3], double acceleration[3])
{
    double rc, p, h, vc;
    int i;
    rc = sqrt(pow(coord[0], 2) + pow(coord[1], 2) + pow(coord[2], 2));

    h = rc - R0 * (1 - A0*pow(coord[2]/rc, 2));
    p = 2e-13 * exp(-((h-200)/60));

    vc = sqrt(pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2));

    for (i = 0; i < 3; i++)
        acceleration[i] = -Sb * p * vc * vel[i];

    return;
}


void get_acceleration_by_solar_pressurel(long double utc_in_mjd, long double coord[3], long double acceleration[3])
{
    long double ro, sun_coord[3];

    long double tdb = tt_to_tdbl(mjd_to_ttl(utc_in_mjd));

    int i;
    get_sun_celestial_positionl(tdb - 0.0057755L, sun_coord);

    ro = sqrtl(powl((sun_coord[0] - coord[0]), 2) + powl((sun_coord[1] - coord[1]), 2) + powl((sun_coord[2] - coord[2]), 2));

    for (i = 0; i < 3; i++)
        acceleration[i] = CREFL * powl(AU/ro, 2) * (coord[i] - sun_coord[i])/ro;

    return;
}


void get_acceleration_by_solar_pressure(double utc_in_mjd, double coord[3], double acceleration[3])
{
    double ro, sun_coord[3];

    double tdb = tt_to_tdb(mjd_to_tt(utc_in_mjd));

    int i;
    get_sun_celestial_position(tdb - 0.0057755, sun_coord);

    ro = sqrt(pow((sun_coord[0] - coord[0]), 2)
              + pow((sun_coord[1] - coord[1]), 2)
              + pow((sun_coord[2] - coord[2]), 2));

    for (i = 0; i < 3; i++)
        acceleration[i] = CREFL * pow(AU/ro, 2) * (coord[i] - sun_coord[i])/ro;

    return;
}


void get_forcesl(long double t, long double coord[3], long double vel[3], long double acceleration[3])
{
    long double acc[3];

    acceleration[0] = 0; acceleration[1] = 0; acceleration[2] = 0;

    get_acceleration_by_earthl(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_moonl(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_sunl(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_atmospherel(coord, vel, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_solar_pressurel(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    return;
}


void get_forces(double t, double coord[3], double vel[3], double acceleration[3])
{
    double acc[3];

    acceleration[0] = 0; acceleration[1] = 0; acceleration[2] = 0;

    get_acceleration_by_earth(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_moon(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_sun(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_atmosphere(coord, vel, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    get_acceleration_by_solar_pressure(t, coord, acc);
    acceleration[0] += acc[0];
    acceleration[1] += acc[1];
    acceleration[2] += acc[2];

    return;
}