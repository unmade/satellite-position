//
// Created by user on 04.06.2015.
//

#include <assert.h>
#include <math.h>
#include "forces_test.h"
#include "../forces/forces.h"

void get_acceleration_by_earthl_test(void)
{
    long double utc_in_mjd = 59865.00000000L;
    long double coord[3], acceleration[3];

    coord[0] = -43203.60633L;
    coord[1] = 932.853125L;
    coord[2] = 105.030427L;

    get_acceleration_by_earthl(utc_in_mjd, coord, acceleration);

    assert(fabsl(fabsl(acceleration[0]) - 0.00021340560717251052L) < 1e-18);
    assert(fabsl(fabsl(acceleration[1]) - 4.6078667089531894e-006L) < 1e-18);
    assert(fabsl(fabsl(acceleration[2]) - 5.1880435966992614e-007L) < 1e-18);

    return;
}

void get_acceleration_by_earth_test(void)
{
    double utc_in_mjd = 59865.00000000;
    double coord[3], acceleration[3];

    coord[0] = -43203.60633;
    coord[1] = 932.853125;
    coord[2] = 105.030427;

    get_acceleration_by_earth(utc_in_mjd, coord, acceleration);

    assert(fabs(fabs(acceleration[0]) - 0.00021340560717251038) < 1e-18);
    assert(fabs(fabs(acceleration[1]) - 4.6078667089531775e-006) < 1e-18);
    assert(fabs(fabs(acceleration[2]) - 5.1880435966992561e-007) < 1e-18);

    return;
}


void get_acceleration_by_moonl_test(void)
{
    long double utc_in_mjd = 59865.0L;
    long double coord[3], acceleration[3];

    coord[0] = -43203.60633L;
    coord[1] = 932.853125L;
    coord[2] = 105.030427L;

    get_acceleration_by_moonl(utc_in_mjd, coord, acceleration);

    assert(fabsl(fabsl(acceleration[0]) - 5.7471499870040367575e-11L) < 1e-18);
    assert(fabsl(fabsl(acceleration[1]) - 4.0468814090081801167e-09L) < 1e-18);
    assert(fabsl(fabsl(acceleration[2]) - 1.8596276722038470893e-09L) < 1e-18);

    return;
}

void get_acceleration_by_moon_test(void)
{
    double utc_in_mjd = 59865.0;
    double coord[3], acceleration[3];

    coord[0] = -43203.60633;
    coord[1] = 932.853125;
    coord[2] = 105.030427;

    get_acceleration_by_moon(utc_in_mjd, coord, acceleration);
    assert(fabs(fabs(acceleration[0]) - 5.7471499863764875e-11) < 1e-18);
    assert(fabs(fabs(acceleration[1]) - 4.0468814090103113e-09) < 1e-18);
    assert(fabs(fabs(acceleration[2]) - 1.8596276722046323e-09) < 1e-18);

    return;
}


void get_acceleration_by_sunl_test(void)
{
    long double utc_in_mjd = 59865.0L;
    long double coord[3], acceleration[3];

    coord[0] = -43203.606331156676L;
    coord[1] = 932.8531254505084L;
    coord[2] = 105.03042685970608L;

    get_acceleration_by_sunl(utc_in_mjd, coord, acceleration);

    assert(fabsl(fabsl(acceleration[0]) - 2.8480208700617356e-009L) < 1e-18);
    assert(fabsl(fabsl(acceleration[1]) - 1.5087331078186935e-009L) < 1e-18);
    assert(fabsl(fabsl(acceleration[2]) - 6.4204449389352687e-010L) < 1e-18);

    return;
}


void get_acceleration_by_sun_test(void)
{
    double utc_in_mjd = 59865.0;
    double coord[3], acceleration[3];

    coord[0] = -43203.606331156676;
    coord[1] = 932.8531254505084;
    coord[2] = 105.03042685970608;

    get_acceleration_by_sun(utc_in_mjd, coord, acceleration);

    assert(fabs(fabs(acceleration[0]) - 2.8480208700633192e-009) < 1e-18);
    assert(fabs(fabs(acceleration[1]) - 1.5087331078190796e-009) < 1e-18);
    assert(fabs(fabs(acceleration[2]) - 6.4204449389369438e-010) < 1e-18);

    return;
}


void test_forces(void)
{
    get_acceleration_by_earthl_test();
    get_acceleration_by_earth_test();

    get_acceleration_by_moonl_test();
    get_acceleration_by_moon_test();

    get_acceleration_by_sunl_test();
    get_acceleration_by_sun_test();

    return;
}