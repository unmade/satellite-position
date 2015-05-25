//
// Created by Леша on 23.05.15.
//

#include <math.h>
#include <stdio.h>
#include "date_converters.h"
#include "../constants.h"
#include "../nutation.h"

#define TT_TAI_DIFF 32.184

double utc_to_mjd(int nyear, int nmonth, int nday,
                  int nhour, int nminute, double nsecond)
{
    int a = (14 - nmonth) / 12;
    int y = nyear + 4800 - a;
    int m = nmonth + 12*a - 3;

    int jdn = nday + ((153*m + 2)/5) + 365*y + (y/4) - (y/100) + (y/400) - 32045;

    double jd = jdn + (((double)nhour-12)/24) + ((double)nminute/1440) + (nsecond/86400);

    return jd - 2400000.5;
}


void mjd_to_utc(double mjd, int *nyear, int *nmonth, int *nday,
                int *nhour, int *nminute, double *nsecond)
{
    double A, B, C, D, E, G, F, I;
    double days, month, year;

    double jd = mjd + 2400000.5 + 0.5;

    F = modf(jd, &I);
    A = trunc((I - 1867216.25) / 36524.25);

    B = (I > 2299160) ? I + 1 + A - trunc(A / 4) : I;
    C = B + 1524;
    D = trunc((C - 122.1) / 365.25);
    E = trunc(365.25 * D);
    G = trunc((C - E) / 30.6001);

    days = C - E + F - trunc(30.6001 * G);
    month = (G < 13) ? G - 1 : G - 13;
    year = (month > 2.5) ? D - 4716 : D - 4715;

    *nyear = (int) year;
    *nmonth = (int) month;
    *nday = (int) days;

    days_to_hms(days - (int)days, nhour, nminute, nsecond);

    return;
}


double utc_to_tt(int nyear, int nmonth, int nday,
                 int nhour, int nminute, double nsecond)
{
    double mjd = utc_to_mjd(nyear, nmonth, nday, nhour, nminute, nsecond);
    double deltaT = get_deltaT(nyear, nmonth);

    return mjd + (deltaT / 86400);
}


double mjd_to_tt(double mjd_in_utc)
{
    int year, month, day, hour, minute;
    double seconds;
    mjd_to_utc(mjd_in_utc, &year, &month, &day, &hour, &minute, &seconds);
    double deltaT = get_deltaT(year, month);

    return mjd_in_utc + (deltaT / 86400);
}


double tt_to_tdb(double tt)
{
    double d = (tt - 51544.5) / 36525;
    double g = (M_PI / 180) * (357.258 + 35999.050 * d);
    return tt + ((0.00168*sin(g+0.0167*sin(g)))/86400);

//    double mjd = utc_to_mjd(2003, 3, 7, 2, 45, 0.0);
//    double t = (mjd + 2400000.5 - 2451545.0) / 36525.0;
//    double dtr = M_PI / 180;
//    double corr;
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


double utc_to_gmst(double utc_in_mjd, double delta_ut)
{
    double ts = utc_in_mjd + delta_ut / 86400;
    double ts_trunc = trunc(ts);
    double s = ts - ts_trunc;
    double tu = (trunc(ts) - 51544.5) / 36525;
    double s0 = 1.753368559233266 + (628.3319706888409
                + (6.770713944903336e-06 - 4.508767234318685e-10*tu)*tu)*tu;
    double freq = 1.002737909350795
                  + (5.900575455674703e-11 - 5.893984333409384e-15*tu)*tu;
    s0 += freq*s*PI2;
    double r = s0 / (PI2);
    double i = trunc(r);

    double gmst = s0 - PI2*i;
    if (gmst < 0)
    {
        gmst += PI2;
    }

    return gmst;
}


double utc_to_gast(int year, int month, int day, int hour,
                   int minute, double second, double delta_ut)
{
    double mjd = utc_to_mjd(year, month, day, hour, minute, second);
    double tt = utc_to_tt(year, month, day, hour, minute, second);
    double tdb = tt_to_tdb(tt);

    double eps = get_eps_mean(tdb);

    double delta_psi, delta_eps;
    get_nutation_parameters(tdb, &delta_psi, &delta_eps);

    double gmst = utc_to_gmst(mjd, delta_ut);

    return gmst + delta_psi*cos(eps);
}


double mjd_to_gast(double utc_in_mjd, double delta_ut)
{
    double tt = mjd_to_tt(utc_in_mjd);
    double tdb = tt_to_tdb(tt);

    double eps = get_eps_mean(tdb);

    double delta_psi, delta_eps;
    get_nutation_parameters(tdb, &delta_psi, &delta_eps);

    double gmst = utc_to_gmst(utc_in_mjd, delta_ut);

    return gmst + delta_psi*cos(eps);
}


void days_to_hms(double days, int *nhour, int *nminute, double *nsecond)
{
    double hour, minute;
    double hours, minutes, seconds;

    hours = days * 24;
    hours = modf(hours, &hour);

    minutes = hours * 60;
    minutes = modf(minutes, &minute);

    seconds = minutes * 60;

    *nhour = (int)hour;
    *nminute = (int)minute;
    *nsecond = seconds;

    return;
}

double get_deltaT(int nyear, int nmonth)
{
    int N = 0;
    switch(nyear)
    {
        case 1972:
            N = (0 < nmonth && nmonth < 7) ? 10 : 11;
            break;
        case 1973:
            N = 12;
            break;
        case 1974:
            N = 13;
            break;
        case 1975:
            N = 14;
            break;
        case 1976:
            N = 15;
            break;
        case 1977:
            N = 16;
            break;
        case 1978:
            N = 17;
            break;
        case 1979:
            N = 18;
            break;
        case 1980:
            N = 19;
            break;
        case 1981:
            N = (0 < nmonth && nmonth < 7) ? 19 : 20;
            break;
        case 1982:
            N = (0 < nmonth && nmonth < 7) ? 20 : 21;
            break;
        case 1983:
            N = (0 < nmonth && nmonth < 7) ? 21 : 22;
            break;
        case 1984:
            N = 22;
            break;
        case 1985:
            N = (0 < nmonth && nmonth < 7) ? 22 : 23;
        case 1986:
        case 1987:
            N = 23;
            break;
        case 1988:
        case 1989:
            N = 24;
            break;
        case 1990:
            N = 25;
            break;
        case 1991:
            N = 26;
            break;
        case 1992:
            N = (0 < nmonth && nmonth < 7) ? 26 : 27;
            break;
        case 1993:
            N = (0 < nmonth && nmonth < 7) ? 27 : 28;
            break;
        case 1994:
            N = (0 < nmonth && nmonth < 7) ? 28 : 29;
            break;
        case 1995:
            N = 29;
            break;
        case 1996:
            N = 30;
        case 1997:
            N = (0 < nmonth && nmonth < 7) ? 30 : 31;
            break;
        case 1998:
            N = 31;
            break;
        case 1999:
        case 2000:
        case 2001:
        case 2002:
        case 2003:
        case 2004:
        case 2005:
            N = 32;
            break;
        case 2006:
        case 2007:
        case 2008:
            N = 33;
            break;
        case 2009:
        case 2010:
        case 2011:
            N = 34;
            break;
        case 2012:
            N = (0 < nmonth && nmonth < 7) ? 34 : 35;
            break;
        case 2013:
        case 2014:
            N = 35;
            break;
        case 2015:
            N = (0 < nmonth && nmonth < 7) ? 35 : 36;
            break;
        default:
            N = 0;
    }

    return N + TT_TAI_DIFF;
}