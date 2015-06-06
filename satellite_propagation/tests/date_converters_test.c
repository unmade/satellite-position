//
// Created by user on 03.06.2015.
//

#include <assert.h>
#include "date_converters_test.h"
#include "../date_converters.h"
#include "../constants.h"

void utc_to_mjdl_test(void)
{
    long double mjd;

    mjd = utc_to_mjdl(2000, 1, 1, 12, 0, 0.0L);
    assert(mjd - MJD2000 < 1e-19);

    mjd = utc_to_mjdl(2036, 2, 29, 2, 45, 0.0L);
    assert((mjd - 64752.114583333336L) < 1e-19);

    return;
}


void utc_to_mjd_test(void)
{
    double mjd;

    mjd = utc_to_mjd(2000, 1, 1, 12, 0, 0.0);
    assert(mjd - MJD2000 < 1e-15);

    mjd = utc_to_mjd(2036, 2, 29, 2, 45, 0.0);
    assert(mjd - 64752.114583333489000 < 1.e-15);

    return;
}


void mjd_to_utcl_test(void)
{
    long double mjd = MJD2000;
    int year, month, day, hour, minute;
    long double seconds;

    mjd_to_utcl(mjd, &year, &month, &day, &hour, &minute, &seconds);
    assert(year == 2000);
    assert(month == 1);
    assert(day == 1);
    assert(hour == 12);
    assert(minute == 0);
    assert(seconds == 0.0L);

    mjd = 64752.114583333336L;
    mjd_to_utcl(mjd, &year, &month, &day, &hour, &minute, &seconds);
    assert(year == 2036);
    assert(month == 2);
    assert(day == 29);
    assert(hour == 2);
    assert(minute == 45);
    assert(seconds < 1e-6L);

    return;
}


void mjd_to_utc_test(void)
{
    double mjd = MJD2000;
    int year, month, day, hour, minute;
    double seconds;

    mjd_to_utc(mjd, &year, &month, &day, &hour, &minute, &seconds);
    assert(year == 2000);
    assert(month == 1);
    assert(day == 1);
    assert(hour == 12);
    assert(minute == 0);
    assert(seconds == 0.0);

    mjd = 64752.114583333336;
    mjd_to_utc(mjd, &year, &month, &day, &hour, &minute, &seconds);
    assert(year == 2036);
    assert(month == 2);
    assert(day == 29);
    assert(hour == 2);
    assert(minute == 45);
    assert(seconds < 1e-4);

    return;
}


void utc_to_ttl_test(void)
{
    long double tt = utc_to_ttl(2003, 11, 15, 15, 35, 0.0L);
    assert(tt - 52958.650048425923L < 1e-11);

    return;
}


void utc_to_tt_test(void)
{
    double tt = utc_to_tt(2003, 11, 15, 15, 35, 0.0);
    assert(tt - 52958.650048425923 < 1e-15);

    return;
}


void mjd_to_ttl_test(void)
{
    long double tt = mjd_to_ttl(utc_to_mjdl(2003, 11, 15, 15, 35, 0.0L));
    assert(tt - 52958.650048425923L < 1e-11);

    return;
}


void mjd_to_tt_test(void)
{
    double tt = mjd_to_tt(utc_to_mjd(2003, 11, 15, 15, 35, 0.0));
    assert(tt - 52958.650048425923 < 1e-15);

    return;
}


void tt_to_tdbl_test(void)
{
    long double tdb = tt_to_tdbl(52958.650048425923L);

    assert(tdb - 52958.650048423013L < 1e-11);

    return;
}


void tt_to_tdb_test(void)
{
    double tdb = tt_to_tdb(52958.650048425923);

    assert(tdb -  52958.650048423013 < 1e-15);

    return;
}


void mjd_to_gmstl_test(void)
{
    long double mjd, gmst;
    long double deltaUT = 0.38415L;

    mjd = utc_to_mjdl(1999, 12, 7, 2, 45, 0.0L);
    gmst = mjd_to_gmstl(mjd, deltaUT);

    assert(fabsl(gmst - 2.0366448505391315L) < 1e-15);
    return;
}


void mjd_to_gmst_test(void)
{
    double mjd, gmst;
    double deltaUT = 0.38415;

    mjd = utc_to_mjd(1999, 12, 7, 2, 45, 0.0);
    gmst = mjd_to_gmst(mjd, deltaUT);

    assert(fabsl(gmst - 2.0366448514991138) < 1e-15);
    return;
}


void utc_to_gastl_test(void)
{
    long double gast = utc_to_gastl(1999, 12, 7, 2, 45, 0.0L, 0.38415L);
    assert(fabsl(gast - 2.0365783990212121L) < 1e-15);
    return;
}


void utc_to_gast_test(void)
{
    double gast = utc_to_gast(1999, 12, 7, 2, 45, 0.0, 0.38415);
    assert(fabs(gast - 2.0365783999811944) < 1e-15);
    return;
}


void mjd_to_gastl_test(void)
{
    long double gast = mjd_to_gastl(utc_to_mjdl(1999, 12, 7, 2, 45, 0.0L), 0.38415L);
    assert(fabsl(gast - 2.0365783990212121L) < 1e-15);
    return;
}


void mjd_to_gast_test(void)
{
    double gast = mjd_to_gast(utc_to_mjd(1999, 12, 7, 2, 45, 0.0), 0.38415);
    assert(fabs(gast - 2.0365783999811944) < 1e-15);
    return;
}


void days_to_hmsl_test(void)
{
    int hour, minute;
    long double seconds;

    days_to_hmsl(0.5L, &hour, &minute, &seconds);
    assert(hour == 12.0L);
    assert(minute == 0.0L);
    assert(seconds == 0.0L);

    days_to_hmsl(0.127L, &hour, &minute, &seconds);
    assert(hour == 3.0L);
    assert(minute == 2.0L);
    assert(fabsl(seconds - 52.8L) < 1e-15);

    return;
}


void days_to_hms_test(void)
{
    int hour, minute;
    double seconds;

    days_to_hms(0.5, &hour, &minute, &seconds);
    assert(hour == 12.0);
    assert(minute == 0.0);
    assert(seconds == 0.0);

    days_to_hms(0.127, &hour, &minute, &seconds);
    assert(hour == 3.0);
    assert(minute == 2.0);
    assert(fabsl(seconds - 52.8) < 1e-12);

    return;
}


void get_deltaTl_test(void)
{
    assert((get_deltaTl(1999, 1)) - 64.184L < 1e-19);
    assert(get_deltaTl(1997, 1) - 63.184L < 1e-19);
    assert(get_deltaTl(1996, 1) - 62.184L < 1e-19);
    return;
}


void test_date_converters(void)
{
    utc_to_mjdl_test();
    utc_to_mjd_test();

    mjd_to_utcl_test();
    mjd_to_utc_test();

    utc_to_ttl_test();
    utc_to_tt_test();
    mjd_to_ttl_test();
    mjd_to_tt_test();

    tt_to_tdbl_test();
    tt_to_tdb_test();

    mjd_to_gmstl_test();
    mjd_to_gmst_test();

    utc_to_gastl_test();
    utc_to_gast_test();
    mjd_to_gastl_test();
    mjd_to_gast_test();


    days_to_hmsl_test();
    days_to_hms_test();

    get_deltaTl_test();
}


//
//
//void days_to_hmsl(long double days, int *nhour, int *nminute, long double *nsecond);
//
//
//void days_to_hms(double days, int *nhour, int *nminute, double *nsecond);
//
//
//long double get_deltaTl(int nyear, int nmonth);
