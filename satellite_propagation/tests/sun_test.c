//
// Created by user on 04.06.2015.
//

#include <assert.h>
#include <math.h>
#include "sun_test.h"
#include "../sun.h"
#include "../constants.h"

void get_sun_ecliptic_positionl_test(void)
{
    long double l, b, r;
    get_sun_ecliptic_positionl(58651.11532621L, &l, &b, &r);

    assert(fabsl(fabsl(l) - 1.4952254608931363L) < 1e-15);
    assert(fabsl(fabsl(b) - 3.5932090028749775e-007L) < 1e-7);
    assert(fabsl(fabsl(r*AU) - 151974058.51469308L) < 1e-7);

    return;
}


void get_sun_ecliptic_position_test(void)
{
    double l, b, r;
    get_sun_ecliptic_position(58651.11532621, &l, &b, &r);
    assert(fabs(fabs(l) - 1.4952254608931084) < 1e-15);
    assert(fabs(fabs(b) - 3.5932090028843685e-007) < 1e-15);
    assert(fabs(fabs(r*AU) - 151974058.51469302) < 1e-7);

    return;
}


void get_sun_celestial_positionl_test(void)
{
    long double coord[3];
    get_sun_celestial_positionl(58651.11532621L, coord);

    coord[0] *= AU;
    coord[1] *= AU;
    coord[2] *= AU;

    assert(fabsl(fabsl(coord[0]) - 12192677.054444425L) < 1e-8);
    assert(fabsl(fabsl(coord[1]) - 138986664.55175757L) < 1e-8);
    assert(fabsl(fabsl(coord[2]) - 60250810.487914562L) < 1e-8);

    return;
}

void get_sun_celestial_position_test(void)
{
    double coord[3];
    get_sun_celestial_position(58651.11532621, coord);

    coord[0] *= AU;
    coord[1] *= AU;
    coord[2] *= AU;

    assert(fabs(fabs(coord[0]) - 12192677.054448647) < 1e-8);
    assert(fabs(fabs(coord[1]) - 138986664.55175722) < 1e-8);
    assert(fabs(fabs(coord[2]) - 60250810.487914428) < 1e-8);

    return;
}


void test_sun(void)
{
    get_sun_ecliptic_positionl_test();
    get_sun_ecliptic_position_test();

    get_sun_celestial_positionl_test();
    get_sun_celestial_position_test();

    return;
}