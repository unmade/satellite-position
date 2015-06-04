//
// Created by user on 04.06.2015.
//

#include <math.h>
#include <assert.h>
#include "moon_test.h"
#include "../moon.h"
#include "../date_converters/date_converters.h"

void get_moon_ecliptic_positionl_test(void)
{
    long double l, b, r;
    get_moon_ecliptic_positionl(58651.11532621L, &l, &b, &r);

    assert(fabsl(fabsl(l) -  4.586800812L) < 1e-8);
    assert(fabsl(fabsl(b) -  0.038541132329576865L) < 1e-17);
    assert(fabsl(fabsl(r) -  388705.81849596079L) < 1e-10);

    return;
}

void get_moon_ecliptic_position_test(void)
{
    double l, b, r;
    get_moon_ecliptic_position(58651.11532621, &l, &b, &r);

    assert(fabs(fabs(l) -  4.5868008123610196) < 1e-15);
    assert(fabs(fabs(b) -  0.038541132329608589) < 1e-17);
    assert(fabs(fabs(r) -  388705.81849595782) < 1e-10);

    return;
}


void get_moon_celestial_positionl_test(void)
{
    long double coord[3];
    get_moon_celestial_positionl(tt_to_tdbl(mjd_to_ttl(utc_to_mjdl(2019, 6, 17, 2, 45, 0.0))), coord);

    assert(fabsl(fabsl(coord[0]) - 50476.051201573231L) < 1e-11);
    assert(fabsl(fabsl(coord[1]) - 359307.79674016603L) < 1e-10);
    assert(fabsl(fabsl(coord[2]) - 139436.08229457741L) < 1e-10);

    return;
}


void get_moon_celestial_position_test(void)
{
    double coord[3];
    get_moon_celestial_position(tt_to_tdb(mjd_to_tt(utc_to_mjd(2019, 6, 17, 2, 45, 0.0))), coord);

    assert(fabs(fabs(coord[0]) - 50476.051187861813) < 1e-11);
    assert(fabs(fabs(coord[1]) - 359307.79674194666) < 1e-10);
    assert(fabs(fabs(coord[2]) - 139436.08229658601) < 1e-10);

    return;
}


void test_moon(void)
{
    get_moon_ecliptic_positionl_test();
    get_moon_ecliptic_position_test();

    get_moon_celestial_positionl_test();
    get_moon_celestial_position_test();

    return;
}