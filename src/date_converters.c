//
// Created by Леша on 23.05.15.
//

#include <math.h>
#include "date_converters.h"
#include "constants.h"
#include "nutation.h"

#define TT_TAI_DIFF 32.18400000000000000000000

long double utc_to_mjdl(int nyear, int nmonth, int nday,
                        int nhour, int nminute, long double nsecond)
{
    int a = (14 - nmonth) / 12;
    int y = nyear + 4800 - a;
    int m = nmonth + 12*a - 3;

    int jdn = nday + ((153*m + 2)/5) + 365*y + (y/4) - (y/100) + (y/400) - 32045;

    long double jd = jdn + ((nhour-12.0L)/24.0L) + (nminute/1440.0L) + (nsecond/86400.0L);

    return jd - DIFF_EPOCH;
}


double utc_to_mjd(int nyear, int nmonth, int nday,
                       int nhour, int nminute, double nsecond)
{
    int a = (14 - nmonth) / 12;
    int y = nyear + 4800 - a;
    int m = nmonth + 12*a - 3;

    int jdn = nday + ((153*m + 2)/5) + 365*y + (y/4) - (y/100) + (y/400) - 32045;

    double jd = jdn + ((nhour-12.0)/24.0) + (nminute/1440.0) + (nsecond/86400.0);

    return jd - DIFF_EPOCH;
}


void mjd_to_utcl(long double mjd, int *nyear, int *nmonth, int *nday,
                 int *nhour, int *nminute, long double *nsecond)
{
    long double A, B, C, D, E, G, F, I;
    long double days, month, year;

    long double jd = mjd + 2400001.0L;

    F = modfl(jd, &I);
    A = truncl((I - 1867216.25L) / 36524.25L);

    B = (I > 2299160.0L) ? I + 1 + A - truncl(A / 4.0L) : I;
    C = B + 1524.0L;
    D = truncl((C - 122.1L) / 365.25L);
    E = truncl(365.25L * D);
    G = truncl((C - E) / 30.6001L);

    days = C - E + F - truncl(30.6001L * G);
    month = (G < 13.0L) ? G - 1 : G - 13.0L;
    year = (month > 2.5L) ? D - 4716.0L : D - 4715.0L;

    *nyear = (int) year;
    *nmonth = (int) month;
    *nday = (int) days;

    days_to_hmsl(days - (int) days, nhour, nminute, nsecond);

    return;
}


void mjd_to_utc(double mjd, int *nyear, int *nmonth, int *nday,
                 int *nhour, int *nminute, double *nsecond)
{
    double A, B, C, D, E, G, F, I;
    double days, month, year;

    double jd = mjd + 2400001.0;

    F = modf(jd, &I);
    A = trunc((I - 1867216.25) / 36524.25);

    B = (I > 2299160) ? I + 1 + A - trunc(A / 4) : I;
    C = B + 1524.0;
    D = trunc((C - 122.1) / 365.25);
    E = trunc(365.25 * D);
    G = trunc((C - E) / 30.6001);

    days = C - E + F - trunc(30.6001 * G);
    month = (G < 13.0) ? G - 1 : G - 13.0;
    year = (month > 2.5) ? D - 4716.0 : D - 4715.0;

    *nyear = (int) year;
    *nmonth = (int) month;
    *nday = (int) days;

    days_to_hms(days - (int) days, nhour, nminute, nsecond);

    return;
}


long double mjd_to_ttl(long double mjd_in_utc)
{
    int year, month, day, hour, minute;
    long double seconds, deltaT;

    mjd_to_utcl(mjd_in_utc, &year, &month, &day, &hour, &minute, &seconds);
    deltaT = get_deltaTl(year, month);

    return mjd_in_utc + (deltaT / SEC_IN_DAY);
}


double mjd_to_tt(double mjd_in_utc)
{
    int year, month, day, hour, minute;
    double seconds, deltaT;

    mjd_to_utc(mjd_in_utc, &year, &month, &day, &hour, &minute, &seconds);
    deltaT = (double)get_deltaTl(year, month);

    return mjd_in_utc + (deltaT / SEC_IN_DAY);
}


long double utc_to_ttl(int nyear, int nmonth, int nday,
                       int nhour, int nminute, long double nsecond)
{
    return mjd_to_ttl(utc_to_mjdl(nyear, nmonth, nday, nhour, nminute, nsecond));
}


double utc_to_tt(int nyear, int nmonth, int nday,
                 int nhour, int nminute, double nsecond)
{
    return mjd_to_tt(utc_to_mjd(nyear, nmonth, nday, nhour, nminute, nsecond));
}


long double tt_to_tdbl(long double tt)
{
    long double d = (tt - MJD2000) / JULIAN_C;
    long double g = GRAD_IN_RAD * (357.258 + 35999.050 * d);
    return tt + ((0.00168 * sinl(g + 0.0167*sinl(g))) / 86400.0);

//    long double mjd = utc_to_mjdl(2003, 3, 7, 2, 45, 0.0);
//    long double t = (mjd + 2400000.5 - 2451545.0) / 36525.0;
//    long double dtr = M_PI / 180;
//    long double corr;
//    corr = 1656.675     * sin(dtr * (35999.3729 * t + 357.5287))
//           + 22.418     * sin(dtr * (32964.467  * t + 246.199))
//           + 13.84      * sin(dtr * (71998.746  * t + 355.057))
//           +  4.77      * sin(dtr * ( 3034.906  * t +  25.463))
//           +  4.677     * sin(dtr * (34777.259  * t + 230.394))
//           + 10.216 * t * sin(dtr * (35999.373  * t + 243.451))
//           +  0.171 * t * sin(dtr * (71998.746  * t + 240.98 ))
//           +  0.027 * t * sin(dtr * ( 1222.114  * t + 194.661))
//           +  0.027 * t * sin(dtr * ( 3034.906  * t + 336.061))
//           +  0.026 * t * sin(dtr * (  -20.186  * t +   9.382))
//           +  0.007 * t * sin(dtr * (29929.562  * t + 264.911))
//           +  0.006 * t * sin(dtr * (  150.678  * t +  59.775))
//           +  0.005 * t * sin(dtr * ( 9037.513  * t + 256.025))
//           +  0.043 * t * sin(dtr * (35999.373  * t + 151.121));
//
//    return tt + (0.000001 * corr / 86400.0);
}


double tt_to_tdb(double tt)
{
    double d = (tt - MJD2000) / JULIAN_C;
    double g = (GRAD_IN_RAD) * (357.258 + 35999.050 * d);

    return tt + (0.00168 * sin(g + 0.0167*sin(g))) / 86400.0;
}


long double mjd_to_gmstl(long double utc_in_mjd, long double delta_ut)
{
    long double ts = utc_in_mjd + delta_ut / SEC_IN_DAY;
    long double ts_trunc = truncl(ts);
    long double s = ts - ts_trunc;
    long double tu = (truncl(ts) - MJD2000) / JULIAN_C;
    long double s0 = 1.753368559233266L + (628.3319706888409L
                                        + (6.770713944903336e-06L - 4.508767234318685e-10L*tu)*tu)*tu;
    long double freq = 1.002737909350795L + (5.900575455674703e-11L - 5.893984333409384e-15L*tu)*tu;

    s0 += freq*s*PI2;
    long double r = s0 / PI2;
    long double i = truncl(r);

    long double gmst = s0 - PI2*i;
    if (gmst < 0)
    {
        gmst += PI2;
    }

    return gmst;
}


double mjd_to_gmst(double utc_in_mjd, double delta_ut)
{
    double ts = utc_in_mjd + delta_ut / SEC_IN_DAY;
    double ts_trunc = trunc(ts);
    double s = ts - ts_trunc;
    double tu = (trunc(ts) - MJD2000) / JULIAN_C;
    double s0 = 1.753368559233266 + (628.3319706888409
                                  + (6.770713944903336e-06 - 4.508767234318685e-10*tu)*tu)*tu;
    double freq = 1.002737909350795 + (5.900575455674703e-11 - 5.893984333409384e-15*tu)*tu;

    s0 += freq*s*PI2;
    double r = s0 / PI2;
    double i = trunc(r);

    double gmst = s0 - PI2*i;
    if (gmst < 0)
    {
        gmst += PI2;
    }

    return gmst;
}


long double mjd_to_gastl(long double utc_in_mjd, long double delta_ut)
{
    long double tt = mjd_to_ttl(utc_in_mjd);
    long double tdb = tt_to_tdbl(tt);

    long double eps = get_eps_meanl(tdb);

    long double delta_psi, delta_eps;
    get_nutation_parametersl(tdb, &delta_psi, &delta_eps);

    long double gmst = mjd_to_gmstl(utc_in_mjd, delta_ut);

    return gmst + delta_psi*cosl(eps);
}


double mjd_to_gast(double utc_in_mjd, double delta_ut)
{
    double tt = mjd_to_tt(utc_in_mjd);
    double tdb = tt_to_tdb(tt);

    double eps = get_eps_mean(tdb);

    double delta_psi, delta_eps;
    get_nutation_parameters(tdb, &delta_psi, &delta_eps);

    double gmst = mjd_to_gmst(utc_in_mjd, delta_ut);

    return gmst + delta_psi*cos(eps);
}


long double utc_to_gastl(int year, int month, int day, int hour,
                         int minute, long double second, long double delta_ut)
{
    return mjd_to_gastl(utc_to_mjdl(year, month, day, hour, minute, second), delta_ut);
}


double utc_to_gast(int year, int month, int day, int hour,
                   int minute, double second, double delta_ut)
{
    return mjd_to_gast(utc_to_mjd(year, month, day, hour, minute, second), delta_ut);
}


void days_to_hmsl(long double days, int *nhour, int *nminute, long double *nsecond)
{
    long double hour, minute;
    long double hours, minutes, seconds;

    hours = days * 24.0L;
    hours = modfl(hours, &hour);

    minutes = hours * 60.0L;
    minutes = modfl(minutes, &minute);

    seconds = minutes * 60.0L;

    *nhour = (int)hour;
    *nminute = (int)minute;
    *nsecond = seconds;

    return;
}


void days_to_hms(double days, int *nhour, int *nminute, double *nsecond)
{
    double hour, minute;
    double hours, minutes, seconds;

    hours = days * 24.0;
    hours = modf(hours, &hour);

    minutes = hours * 60.0;
    minutes = modf(minutes, &minute);

    seconds = minutes * 60.0;

    *nhour = (int)hour;
    *nminute = (int)minute;
    *nsecond = seconds;

    return;
}


long double get_deltaTl(int nyear, int nmonth)
{
    int leap_second = 0;

    switch(nyear)
    {
        case 1972:
            leap_second = (0 < nmonth && nmonth < 7) ? 10 : 11;
            break;
        case 1973:
            leap_second = 12;
            break;
        case 1974:
            leap_second = 13;
            break;
        case 1975:
            leap_second = 14;
            break;
        case 1976:
            leap_second = 15;
            break;
        case 1977:
            leap_second = 16;
            break;
        case 1978:
            leap_second = 17;
            break;
        case 1979:
            leap_second = 18;
            break;
        case 1980:
            leap_second = 19;
            break;
        case 1981:
            leap_second = (0 < nmonth && nmonth < 7) ? 19 : 20;
            break;
        case 1982:
            leap_second = (0 < nmonth && nmonth < 7) ? 20 : 21;
            break;
        case 1983:
            leap_second = (0 < nmonth && nmonth < 7) ? 21 : 22;
            break;
        case 1984:
            leap_second = 22;
            break;
        case 1985:
            leap_second = (0 < nmonth && nmonth < 7) ? 22 : 23;
            break;
        case 1986:
        case 1987:
            leap_second = 23;
            break;
        case 1988:
        case 1989:
            leap_second = 24;
            break;
        case 1990:
            leap_second = 25;
            break;
        case 1991:
            leap_second = 26;
            break;
        case 1992:
            leap_second = (0 < nmonth && nmonth < 7) ? 26 : 27;
            break;
        case 1993:
            leap_second = (0 < nmonth && nmonth < 7) ? 27 : 28;
            break;
        case 1994:
            leap_second = (0 < nmonth && nmonth < 7) ? 28 : 29;
            break;
        case 1995:
            leap_second = 29;
            break;
        case 1996:
            leap_second = 30;
            break;
        case 1997:
            leap_second = (0 < nmonth && nmonth < 7) ? 30 : 31;
            break;
        case 1998:
            leap_second = 31;
            break;
        case 1999:
        case 2000:
        case 2001:
        case 2002:
        case 2003:
        case 2004:
        case 2005:
            leap_second = 32;
            break;
        case 2006:
        case 2007:
        case 2008:
            leap_second = 33;
            break;
        case 2009:
        case 2010:
        case 2011:
            leap_second = 34;
            break;
        case 2012:
            leap_second = (0 < nmonth && nmonth < 7) ? 34 : 35;
            break;
        case 2013:
        case 2014:
            leap_second = 35;
            break;
        case 2015:
            leap_second = (0 < nmonth && nmonth < 7) ? 35 : 36;
            break;
        default:
            leap_second = (nyear > 1972) ? 36 : 10;
//            leap_second = (nyear > 1972) ? 32 : 10;
    }

    return leap_second + TT_TAI_DIFF;
}