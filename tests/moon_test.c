//
// Created by user on 04.06.2015.
//

#include <math.h>
#include <assert.h>
#include "moon_test.h"
#include "moon.h"
#include "date_converters.h"

void get_moon_ecliptic_positionl_test(void)
{
    long double l, b, r;
    get_moon_ecliptic_positionl(58651.11532621L, &l, &b, &r);

    assert(fabsl(fabsl(l) - 1638.2149806790538514311705853288004L) < 1e-8);
    assert(fabsl(fabsl(b) - 0.038541132329579107003994986255235311L) < 1e-17);
    assert(fabsl(fabsl(r) - 388705.59300521719634957662492524832L) < 1e-10);

    return;
}

void get_moon_ecliptic_position_test(void)
{
    double l, b, r;
    get_moon_ecliptic_position(58651.11532621, &l, &b, &r);

    assert(fabs(fabs(l) - 1638.2149806790535) < 1e-11);
    assert(fabs(fabs(b) - 0.038541132329611073) < 1e-11);
    assert(fabs(fabs(r) - 388705.59300520923) < 1e-7);

    return;
}


void get_moon_celestial_positionl_test(void)
{
    long double coord[3];
    get_moon_celestial_positionl(tt_to_tdbl(mjd_to_ttl(utc_to_mjdl(2019, 6, 17, 2, 45, 0.0))), coord);

    assert(fabsl(fabsl(coord[0]) - 50476.021920101909863376477005658671L) < 1e-15);
    assert(fabsl(fabsl(coord[1]) - 359307.58830340288884030996996443719L) < 1e-15);
    assert(fabsl(fabsl(coord[2]) - 139436.00140680852236130249366397038L) < 1e-15);

    return;
}


void get_moon_celestial_position_test(void)
{
    double coord[3];
    get_moon_celestial_position(tt_to_tdb(mjd_to_tt(utc_to_mjd(2019, 6, 17, 2, 45, 0.0))), coord);

    assert(fabs(fabs(coord[0]) - 50476.021906370261) < 1e-9);
    assert(fabs(fabs(coord[1]) - 359307.58830518054) < 1e-8);
    assert(fabs(fabs(coord[2]) - 139436.0014088208) < 1e-8);

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