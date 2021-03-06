//
// Created by Леша on 24.05.15.
//

#include "nutation.h"
#include "constants.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"


// 1980 IAU Theory of Nutation  ( J. M. Wahr )
static const short NUT_ARGS[106][5] =
        { {   1 ,    0 ,    0 ,    0 ,    0 } ,
          {   2 ,    0 ,    0 ,    0 ,    0 } ,
          {   1 ,   -2 ,    0 ,    2 ,    0 } ,
          {   0 ,    2 ,    0 ,   -2 ,    0 } ,
          {   2 ,   -2 ,    0 ,    2 ,    0 } ,
          {   0 ,    1 ,   -1 ,    0 ,   -1 } ,
          {   1 ,    0 ,   -2 ,    2 ,   -2 } ,
          {   1 ,    2 ,    0 ,   -2 ,    0 } ,
          {   2 ,    0 ,    0 ,    2 ,   -2 } ,
          {   0 ,    0 ,    1 ,    0 ,    0 } ,
          {   2 ,    0 ,    1 ,    2 ,   -2 } ,
          {   2 ,    0 ,   -1 ,    2 ,   -2 } ,
          {   1 ,    0 ,    0 ,    2 ,   -2 } ,
          {   0 ,    2 ,    0 ,    0 ,   -2 } ,
          {   0 ,    0 ,    0 ,    2 ,   -2 } ,
          {   0 ,    0 ,    2 ,    0 ,    0 } ,
          {   1 ,    0 ,    1 ,    0 ,    0 } ,
          {   2 ,    0 ,    2 ,    2 ,   -2 } ,
          {   1 ,    0 ,   -1 ,    0 ,    0 } ,
          {   1 ,   -2 ,    0 ,    0 ,    2 } ,
          {   1 ,    0 ,   -1 ,    2 ,   -2 } ,
          {   1 ,    2 ,    0 ,    0 ,   -2 } ,
          {   1 ,    0 ,    1 ,    2 ,   -2 } ,
          {   0 ,    1 ,    0 ,    0 ,   -1 } ,
          {   0 ,    2 ,    1 ,    0 ,   -2 } ,
          {   1 ,    0 ,    0 ,   -2 ,    2 } ,
          {   0 ,    0 ,    1 ,   -2 ,    2 } ,
          {   2 ,    0 ,    1 ,    0 ,    0 } ,
          {   1 ,   -1 ,    0 ,    0 ,    1 } ,
          {   0 ,    0 ,    1 ,    2 ,   -2 } ,
          {   2 ,    0 ,    0 ,    2 ,    0 } ,
          {   0 ,    1 ,    0 ,    0 ,    0 } ,
          {   1 ,    0 ,    0 ,    2 ,    0 } ,
          {   2 ,    1 ,    0 ,    2 ,    0 } ,
          {   0 ,    1 ,    0 ,    0 ,   -2 } ,
          {   2 ,   -1 ,    0 ,    2 ,    0 } ,
          {   0 ,    0 ,    0 ,    0 ,    2 } ,
          {   1 ,    1 ,    0 ,    0 ,    0 } ,
          {   1 ,   -1 ,    0 ,    0 ,    0 } ,
          {   2 ,   -1 ,    0 ,    2 ,    2 } ,
          {   1 ,    1 ,    0 ,    2 ,    0 } ,
          {   2 ,    0 ,    0 ,    2 ,    2 } ,
          {   0 ,    2 ,    0 ,    0 ,    0 } ,
          {   2 ,    1 ,    0 ,    2 ,   -2 } ,
          {   2 ,    2 ,    0 ,    2 ,    0 } ,
          {   0 ,    0 ,    0 ,    2 ,    0 } ,
          {   1 ,   -1 ,    0 ,    2 ,    0 } ,
          {   1 ,   -1 ,    0 ,    0 ,    2 } ,
          {   1 ,    1 ,    0 ,    0 ,   -2 } ,
          {   1 ,   -1 ,    0 ,    2 ,    2 } ,
          {   0 ,    1 ,    1 ,    0 ,   -2 } ,
          {   2 ,    0 ,    1 ,    2 ,    0 } ,
          {   2 ,    0 ,   -1 ,    2 ,    0 } ,
          {   2 ,    1 ,    0 ,    2 ,    2 } ,
          {   0 ,    1 ,    0 ,    0 ,    2 } ,
          {   2 ,    2 ,    0 ,    2 ,   -2 } ,
          {   1 ,    0 ,    0 ,    0 ,    2 } ,
          {   1 ,    0 ,    0 ,    2 ,    2 } ,
          {   1 ,    1 ,    0 ,    2 ,   -2 } ,
          {   1 ,    0 ,    0 ,    0 ,   -2 } ,
          {   0 ,    1 ,   -1 ,    0 ,    0 } ,
          {   1 ,    2 ,    0 ,    2 ,    0 } ,
          {   0 ,    0 ,    1 ,    0 ,   -2 } ,
          {   0 ,    1 ,    0 ,   -2 ,    0 } ,
          {   0 ,    0 ,    0 ,    0 ,    1 } ,
          {   0 ,    1 ,    1 ,    0 ,    0 } ,
          {   0 ,    1 ,    0 ,    2 ,    0 } ,
          {   2 ,    1 ,   -1 ,    2 ,    0 } ,
          {   2 ,   -1 ,   -1 ,    2 ,    2 } ,
          {   1 ,   -2 ,    0 ,    0 ,    0 } ,
          {   2 ,    3 ,    0 ,    2 ,    0 } ,
          {   2 ,    0 ,   -1 ,    2 ,    2 } ,
          {   2 ,    1 ,    1 ,    2 ,    0 } ,
          {   1 ,   -1 ,    0 ,    2 ,   -2 } ,
          {   1 ,    2 ,    0 ,    0 ,    0 } ,
          {   2 ,    1 ,    0 ,    0 ,    0 } ,
          {   0 ,    3 ,    0 ,    0 ,    0 } ,
          {   2 ,    0 ,    0 ,    2 ,    1 } ,
          {   2 ,   -1 ,    0 ,    0 ,    0 } ,
          {   0 ,    1 ,    0 ,    0 ,   -4 } ,
          {   2 ,   -2 ,    0 ,    2 ,    2 } ,
          {   2 ,   -1 ,    0 ,    2 ,    4 } ,
          {   0 ,    2 ,    0 ,    0 ,   -4 } ,
          {   2 ,    1 ,    1 ,    2 ,   -2 } ,
          {   1 ,    1 ,    0 ,    2 ,    2 } ,
          {   2 ,   -2 ,    0 ,    2 ,    4 } ,
          {   2 ,   -1 ,    0 ,    4 ,    0 } ,
          {   0 ,    1 ,   -1 ,    0 ,   -2 } ,
          {   1 ,    2 ,    0 ,    2 ,   -2 } ,
          {   2 ,    2 ,    0 ,    2 ,    2 } ,
          {   1 ,    1 ,    0 ,    0 ,    2 } ,
          {   2 ,    0 ,    0 ,    4 ,   -2 } ,
          {   2 ,    3 ,    0 ,    2 ,   -2 } ,
          {   0 ,    1 ,    0 ,    2 ,   -2 } ,
          {   1 ,    0 ,    1 ,    2 ,    0 } ,
          {   1 ,   -1 ,   -1 ,    0 ,    2 } ,
          {   1 ,    0 ,    0 ,   -2 ,    0 } ,
          {   2 ,    0 ,    0 ,    2 ,   -1 } ,
          {   0 ,    0 ,    1 ,    0 ,    2 } ,
          {   0 ,    1 ,    0 ,   -2 ,   -2 } ,
          {   1 ,    0 ,   -1 ,    2 ,    0 } ,
          {   1 ,    1 ,    1 ,    0 ,   -2 } ,
          {   0 ,    1 ,    0 ,   -2 ,    2 } ,
          {   0 ,    2 ,    0 ,    0 ,    2 } ,
          {   2 ,    0 ,    0 ,    2 ,    4 } ,
          {   0 ,    0 ,    1 ,    0 ,    1 } } ;

// 1980 IAU Theory of Nutation  ( J. M. Wahr )
// TODO: сделать double?
static const float NUT_COEFF[106][4] =
        { {   -17.199600 ,     -0.017420 ,      9.202500 ,      0.000890 } ,
          {     0.206200 ,      0.000020 ,     -0.089500 ,      0.000050 } ,
          {     0.004600 ,      0.000000 ,     -0.002400 ,      0.000000 } ,
          {     0.001100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000200 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -1.318700 ,     -0.000160 ,      0.573600 ,     -0.000310 } ,
          {     0.142600 ,     -0.000340 ,      0.005400 ,     -0.000010 } ,
          {    -0.051700 ,      0.000120 ,      0.022400 ,     -0.000060 } ,
          {     0.021700 ,     -0.000050 ,     -0.009500 ,      0.000030 } ,
          {     0.012900 ,      0.000010 ,     -0.007000 ,      0.000000 } ,
          {     0.004800 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.002200 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.001700 ,     -0.000010 ,      0.000000 ,      0.000000 } ,
          {    -0.001500 ,      0.000000 ,      0.000900 ,      0.000000 } ,
          {    -0.001600 ,      0.000010 ,      0.000700 ,      0.000000 } ,
          {    -0.001200 ,      0.000000 ,      0.000600 ,      0.000000 } ,
          {    -0.000600 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {    -0.000500 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {     0.000400 ,      0.000000 ,     -0.000200 ,      0.000000 } ,
          {     0.000400 ,      0.000000 ,     -0.000200 ,      0.000000 } ,
          {    -0.000400 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.227400 ,     -0.000020 ,      0.097700 ,     -0.000050 } ,
          {     0.071200 ,      0.000010 ,     -0.000700 ,      0.000000 } ,
          {    -0.038600 ,     -0.000040 ,      0.020000 ,      0.000000 } ,
          {    -0.030100 ,      0.000000 ,      0.012900 ,     -0.000010 } ,
          {    -0.015800 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {     0.012300 ,      0.000000 ,     -0.005300 ,      0.000000 } ,
          {     0.006300 ,      0.000000 ,     -0.000200 ,      0.000000 } ,
          {     0.006300 ,      0.000010 ,     -0.003300 ,      0.000000 } ,
          {    -0.005800 ,     -0.000010 ,      0.003200 ,      0.000000 } ,
          {    -0.005900 ,      0.000000 ,      0.002600 ,      0.000000 } ,
          {    -0.005100 ,      0.000000 ,      0.002700 ,      0.000000 } ,
          {    -0.003800 ,      0.000000 ,      0.001600 ,      0.000000 } ,
          {     0.002900 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {     0.002900 ,      0.000000 ,     -0.001200 ,      0.000000 } ,
          {    -0.003100 ,      0.000000 ,      0.001300 ,      0.000000 } ,
          {     0.002600 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {     0.002100 ,      0.000000 ,     -0.001000 ,      0.000000 } ,
          {     0.001600 ,      0.000000 ,     -0.000800 ,      0.000000 } ,
          {    -0.001300 ,      0.000000 ,      0.000700 ,      0.000000 } ,
          {    -0.001000 ,      0.000000 ,      0.000500 ,      0.000000 } ,
          {    -0.000700 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000700 ,      0.000000 ,     -0.000300 ,      0.000000 } ,
          {    -0.000700 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {    -0.000800 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {     0.000600 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000600 ,      0.000000 ,     -0.000300 ,      0.000000 } ,
          {    -0.000600 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {    -0.000700 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {     0.000600 ,      0.000000 ,     -0.000300 ,      0.000000 } ,
          {    -0.000500 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {     0.000500 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000500 ,      0.000000 ,      0.000300 ,      0.000000 } ,
          {    -0.000400 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000400 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000400 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000300 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000200 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000300 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {     0.000200 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000200 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {     0.000200 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000200 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {     0.000200 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000200 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000200 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000100 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,     -0.000100 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {    -0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } ,
          {     0.000100 ,      0.000000 ,      0.000000 ,      0.000000 } } ;


long double get_eps_meanl(long double tdb)
{
    long double ts = (tdb - MJD2000) / JULIAN_C;
    long double ts2 = ts*ts;
    long double ts3 = ts*ts2;
    return SEC_IN_RAD*(84381.448L - 46.8150L*ts - 0.00059L*ts2 + 0.001813L*ts3);
}

double get_eps_mean(double tdb)
{
    double ts = (tdb - MJD2000) / JULIAN_C;
    double ts2 = ts*ts;
    double ts3 = ts*ts2;
    return SEC_IN_RAD*(84381.448 - 46.8150*ts - 0.00059*ts2 + 0.001813*ts3);
}


void get_fund_argsl(long double tdb, long double fund_args[5])
{
    long double dt = (tdb - MJD2000) / JULIAN_C;
    long double dt2 = dt*dt;
    long double dt3 = dt*dt2;

    fund_args[0] = (218.31643250L + 481267.8812772222L*dt
                    -0.00161167L*dt2 + 0.00000528L*dt3) * GRAD_IN_RAD ;
    fund_args[1] = (134.96298139L + 477198.8673980556L*dt
                    +0.00869722L*dt2 + 0.00001778L*dt3) * GRAD_IN_RAD;
    fund_args[2] = (357.52772333L + 35999.05034L*dt
                    -0.00016028e0L*dt2 - 0.00000333L*dt3) * GRAD_IN_RAD;
    fund_args[3] = (93.27191028L + 483202.0175380555L*dt
                    -0.00368250L*dt2 + 0.00000306L*dt3) * GRAD_IN_RAD;
    fund_args[4] = (297.85036306L + 445267.11148L*dt
                    -0.00191417L*dt2 + 0.00000528L*dt3) * GRAD_IN_RAD;

    return;
}

void get_fund_args(double tdb, double fund_args[5])
{
    double dt = (tdb - MJD2000) / JULIAN_C;
    double dt2 = dt*dt;
    double dt3 = dt*dt2;
    fund_args[0] = (218.31643250 + 481267.8812772222*dt
                   -0.00161167*dt2 + 0.00000528*dt3) * GRAD_IN_RAD;
    fund_args[1] = (134.96298139 + 477198.8673980556*dt
                   +0.00869722*dt2 + 0.00001778*dt3) * GRAD_IN_RAD;
    fund_args[2] = (357.52772333 + 35999.05034*dt
                   -0.00016028e0*dt2 - 0.00000333*dt3) * GRAD_IN_RAD;
    fund_args[3] = (93.27191028 + 483202.0175380555*dt
                   -0.00368250*dt2 + 0.00000306*dt3) * GRAD_IN_RAD;
    fund_args[4] = (297.85036306 + 445267.11148*dt
                   -0.00191417*dt2 + 0.00000528*dt3) * GRAD_IN_RAD;

    return;
}


void get_corr_fund_argsl(long double tdb, long double corr_fund_args[5])
{
    long double s[7];
    long double ts, dr;

    get_fund_argsl(tdb, corr_fund_args);
    
    ts = (tdb - MJD2000) / JULIAN_C;
    s[0] = sinl(1.24614L + 0.35255L*ts);
    s[1] = sinl(1.75106L + 0.28325L*ts);
    s[2] = sinl(1.05727L - 2.31868L*ts);
    s[3] = sinl(2.18240L - 33.7571L*ts);
    s[4] = sinl(0.65961L - 33.79719L*ts);
    s[5] = sinl(2.68173L - 2.62983L*ts);
    s[6] = sinl(0.93890L - 33.77281L*ts);

    dr = 0.31L*s[1] + 14.27L*s[2];
    corr_fund_args[0] += SEC_IN_RAD*(0.84L*s[0] + dr + 7.26L*s[3] + 0.28L*s[4] + 0.24L*s[5]);        // { lambda rd }
    corr_fund_args[1] += SEC_IN_RAD*(2.94L*s[0] + dr + 9.34L*s[3]+1.12L*s[4] + 0.83L*s[5]);          // { l  radian }
    corr_fund_args[2] -= SEC_IN_RAD*(6.4L*s[0] + 1.89L*s[5]);                                        // { l' radian }
    corr_fund_args[3] += SEC_IN_RAD*(0.21L*s[0]+dr-88.7L*s[3]-15.3L*s[4] + 0.24L*s[5] - 1.86L*s[6]); // { F }
    corr_fund_args[4] += SEC_IN_RAD*(7.24L*s[0]+dr + 7.26L*s[3] + 0.28L*s[4]+2.13L*s[5]);             // { D  radian }
}


void get_corr_fund_args(double tdb, double corr_fund_args[5])
{
    double s[7];
    double ts, dr;

    get_fund_args(tdb, corr_fund_args);

    ts = (tdb - MJD2000) / JULIAN_C;
    s[0] = sin(1.24614 + 0.35255*ts);
    s[1] = sin(1.75106 + 0.28325*ts);
    s[2] = sin(1.05727 - 2.31868*ts);
    s[3] = sin(2.18240 - 33.7571*ts);
    s[4] = sin(0.65961 - 33.79719*ts);
    s[5] = sin(2.68173 - 2.62983*ts);
    s[6] = sin(0.93890 - 33.77281*ts);

    dr = 0.31*s[1] + 14.27*s[2];
    corr_fund_args[0] += SEC_IN_RAD*(0.84*s[0] + dr + 7.26*s[3] + 0.28*s[4] + 0.24*s[5]);       // { lambda rd }
    corr_fund_args[1] += SEC_IN_RAD*(2.94*s[0] + dr + 9.34*s[3]+1.12*s[4] + 0.83*s[5]);         // { l  radian }
    corr_fund_args[2] -= SEC_IN_RAD*(6.4*s[0] + 1.89*s[5]);                                     // { l' radian }
    corr_fund_args[3] += SEC_IN_RAD*(0.21*s[0]+dr-88.7*s[3]-15.3*s[4] + 0.24*s[5] - 1.86*s[6]); // { F }
    corr_fund_args[4] += SEC_IN_RAD*(7.24*s[0]+dr + 7.26*s[3] + 0.28*s[4]+2.13*s[5]);           // { D  radian }
}


void get_nutation_parametersl(long double tdb, long double *psi, long double *eps)
{

    long double ts = (tdb - MJD2000) / JULIAN_C;
    long double delta_psi = 0, delta_eps = 0;
    int i, j;

    long double fund_args[5];
    get_fund_argsl(tdb, fund_args);

    long double r = 0;
    for (i = 0; i < 106; i++)
    {
        r = NUT_ARGS[i][0] * (fund_args[0] - fund_args[3]);
        for (j = 1; j < 5; j++)
        {
            r += NUT_ARGS[i][j] * fund_args[j];
        }
        delta_psi += (NUT_COEFF[i][0] + ts* NUT_COEFF[i][1]) * sinl(r);
        delta_eps += (NUT_COEFF[i][2] + ts* NUT_COEFF[i][3]) * cosl(r);
    }
    *psi = delta_psi * SEC_IN_RAD;
    *eps = delta_eps * SEC_IN_RAD;

    return;
}

void get_nutation_parameters(double tdb, double *psi, double *eps)
{
    long double ts = (tdb - MJD2000) / JULIAN_C;
    long double delta_psi = 0, delta_eps = 0;
    int i, j;

    long double fund_args[5];
    get_fund_argsl(tdb, fund_args);

    long double r = 0;
    for (i = 0; i < 106; i++)
    {
        r = NUT_ARGS[i][0] * (fund_args[0] - fund_args[3]);
        for (j = 1; j < 5; j++)
        {
            r += NUT_ARGS[i][j] * fund_args[j];
        }
        delta_psi += (NUT_COEFF[i][0] + ts* NUT_COEFF[i][1]) * sinl(r);
        delta_eps += (NUT_COEFF[i][2] + ts* NUT_COEFF[i][3]) * cosl(r);
    }
    *psi = (double)(delta_psi * SEC_IN_RAD);
    *eps = (double)(delta_eps * SEC_IN_RAD);

    return;
}


void get_nutation_matrixl(long double tdb, long double nutation_matrix[3][3])
{
    long double Rx_deps[3][3], Rz_psi[3][3], Rx_eps[3][3], R[3][3];
    long double delta_psi, delta_eps;
    long double eps = get_eps_meanl(tdb);

    get_nutation_parametersl(tdb, &delta_psi, &delta_eps);

    rotxl(-eps - delta_eps, Rx_deps);
    rotzl(-delta_psi, Rz_psi);
    rotxl(eps, Rx_eps);
    mult_matricesl(Rz_psi, Rx_eps, R);
    mult_matricesl(Rx_deps, R, nutation_matrix);

    return;
}

void get_nutation_matrix(double tdb, double nutation_matrix[3][3])
{
    double Rx_deps[3][3], Rz_psi[3][3], Rx_eps[3][3], R[3][3];
    double delta_psi, delta_eps;
    double eps = get_eps_mean(tdb);

    get_nutation_parameters(tdb, &delta_psi, &delta_eps);

    rotx(-eps - delta_eps, Rx_deps);
    rotz(-delta_psi, Rz_psi);
    rotx(eps, Rx_eps);
    mult_matrices(Rz_psi, Rx_eps, R);
    mult_matrices(Rx_deps, R, nutation_matrix);

    return;
}