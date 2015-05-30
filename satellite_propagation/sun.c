//
// Created by Леша on 29.05.15.
//

#include "sun.h"
#include "constants.h"
#include "coordinates_converters.h"
#include "nutation.h"
#include "rotation_matrix.h"
#include "precession.h"
#include "matrix_operations.h"

long double frac(long double ld)
{
    return ld - (int)ld;
}

void get_sun_ecliptic_position(long double tdb, long double *l, long double *b, long double *r)
{
    long double ts;
    long double um2, um3, um4, um5, um6; // Средние аномалии планет
    long double dal, daf, dad; // Средние аргументы модели движения Луны
    long double dl, dr, db;
    long double u, c, s;
    
    ts = (tdb - 51544.5) / 36525;

    um2 = PI2*frac(0.1387306 + 162.5485917*ts); // { the Venus   M2 }
    um3 = PI2*frac(0.9931266 + 99.9973604*ts);  // { the Earth   M3 }
    um4 = PI2*frac(0.0543250 + 53.1666028*ts);  // { the Mars    M4 }
    um5 = PI2*frac(0.0551750 + 8.4293972*ts);   // { the Jupiter M5 }
    um6 = PI2*frac(0.8816500 + 3.3938722*ts);   // { the Saturn  M6 }

    // { mean arguments of lunar orbit in radian }
    dad = PI2*frac(0.8274+1236.8531*ts); // { elongation of the Moon    D  D }
    dal = PI2*frac(0.3749+1325.5524*ts); // { mean anomaly of the Moon  l  A }
    daf = PI2*frac(0.2591+1342.2278*ts); // { mean argument of latitude F  U }

    //{ keplerian terms and perturbations by Venus }
    dl = 0.0 ; dr = 0.0; db = 0.0; // { initial nullo corrections }
    u = um3;   //{ 1 0 0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.22*c+6892.76*s;
    dr = dr-16707.37*c-0.54*s;
    dl = dl+(-0.06*c-17.35*s)*ts;      //{ 1 0 1 }
    dr = dr+(+42.04*c-0.15*s)*ts;
    dl = dl+(-0.01*c-0.05*s)*powl(ts, 2); //{ 1 0 2 }
    dr = dr+(+0.13*c-0.02*s)*powl(ts, 2);
    u = 2*um3;  //{ 2 0 1 }
    c = cosl(u); s = sinl(u);
    dl = dl+71.98*s;
    dr = dr-139.57*c;
    dl = dl-0.36*s*ts;
    dr = dr+0.70*c*ts;
    u = 3*um3; ///{ 3 0 0 }
    c = cosl(u); s = sinl(u);
    dl = dl+1.04*s;
    dr = dr-1.75*c;
    u = -um2;  //{ 0-1 0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.03*c-0.07*s;
    dr = dr-0.16*c-0.07*s;
    db = db+0.02*c-0.02*s;
    u = um3-um2;// { 1,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+2.35*c-4.23*s;
    dr = dr-4.75*c-2.64*s;
    u = um3-2*um2;// { 1,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.10*c+0.06*s;
    dr = dr+0.12*c+0.20*s;
    db = db+0.02*c;
    u = 2*um3-um2; //{ 2,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.06*c-0.03*s;
    dr = dr+0.20*c-0.01*s;
    db = db+0.01*c-0.09*s;
    u = 2*um3-2*um2;//{ 2,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-4.70*c+2.90*s;
    dr = dr+8.28*c+13.42*s;
    db = db+0.01*c-0.01*s;
    u = 3*um3-2*um2; //{ 3,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+1.80*c-1.74*s;
    dr = dr-1.44*c-1.57*s;
    db = db+0.04*c-0.06*s;
    u = 3*um3-3*um2; //{ 3,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.67*c+0.03*s;
    dr = dr+0.11*c+2.43*s;
    db = db+0.01*c;
    u = 4*um3-2*um2; //{ 4,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.03*c-0.03*s;
    dr = dr+0.10*c+0.09*s;
    db = db+0.01*c-0.01*s;
    u = 4*um3-3*um2; //{ 4,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+1.51*c-0.40*s;
    dr = dr-0.88*c-3.36*s;
    db = db+0.18*c-0.10*s;
    u = 4*um3-4*um2;// { 4,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.19*c-0.09*s;
    dr = dr-0.38*c+0.77*s;
    u = 5*um3-3*um2;// { 5,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.76*c-0.68*s;
    dr = dr+0.30*c+0.37*s;
    db = db+0.01*c;
    u = 5*um3-4*um2; //{ 5,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.14*c-0.04*s;
    dr = dr-0.11*c+0.43*s;
    db = db-0.03*c;
    u = 5*um3-5*um2; //{ 5,-5,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.05*c-0.07*s;
    dr = dr-0.31*c+0.21*s;
    u = 6*um3-4*um2; //{ 6,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.15*c-0.04*s;
    dr = dr-0.06*c-0.21*s;
    db = db+0.01*c;
    u = 6*um3-5*um2; //{ 6,-5,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.03*c-0.03*s;
    dr = dr-0.09*c+0.09*s;
    db = db-0.01*c;
    u = 6*um3-6*um2; //{ 6,-6,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.04*s;
    dr = dr-0.18*c+0.02*s;
    u = 7*um3-5*um2; //{ 7,-5,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.12*c-0.03*s;
    dr = dr-0.08*c+0.31*s;
    db = db-0.02*c-0.01*s;
    
    //{ perturbations by Mars }
    u = um3-um4; //{ 1,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.22*c+0.17*s;
    dr = dr-0.21*c-0.27*s;
    u = um3-2*um4;// { 1,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-1.66*c+0.62*s;
    dr = dr+0.16*c+0.28*s;
    u = 2*um3-2*um4; //{ 2,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+1.96*c+0.57*s;
    dr = dr-1.32*c+4.55*s;
    db = db+0.01*s;
    u = 2*um3-3*um4;// { 2,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.40*c+0.15*s;
    dr = dr-0.17*c+0.46*s;
    u = 2*um3-4*um4; //{ 2,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.53*c+0.26*s;
    dr = dr+0.09*c-0.22*s;
    u = 3*um3-3*um4;// { 3,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.05*c+0.12*s;
    dr = dr-0.35*c+0.15*s;
    u = 3*um3-4*um4;// { 3,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.13*c-0.48*s;
    dr = dr+1.06*c-0.29*s;
    db = db+0.01*c;
    u = 3*um3-5*um4; //{ 3,-5,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.04*c-0.20*s;
    dr = dr+0.20*c-0.04*s;
    u = 4*um3-4*um4;// { 4,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.03*s;
    dr = dr+0.10*c+0.04*s;
    u = 4*um3-5*um4;// { 4,-5,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.05*c-0.07*s;
    dr = dr+0.20*c+0.14*s;
    u = 4*um3-6*um4; //{ 4,-6,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.10*c+0.11*s;
    dr = dr-0.23*c-0.22*s;
    u = 5*um3-7*um4; //{ 5,-7,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.05*c;
    dr = dr+0.01*c-0.14*s;
    u = 5*um3-8*um4; //{ 5,-8,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.05*c+0.01*s;
    dr = dr-0.02*c+0.10*s;
    //{ perturbations by Jupiter }
    u = -um3-um5;// {-1,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.01*c+0.07*s;
    dr = dr+0.18*c-0.02*s;
    db = db-0.02*s;
    u = -um5; //{ 0,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.31*c+2.58*s;
    dr = dr+0.52*c+0.34*s;
    db = db+0.02*c;
    u = um3-um5; //{ 1,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-7.21*c-0.06*s;
    dr = dr+0.13*c-16.27*s;
    db = db-0.02*s;
    u = um3-2*um5; //{ 1,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.54*c-1.52*s;
    dr = dr+3.09*c-1.12*s;
    db = db+0.01*c-0.17*s;
    u = um3-3*um5;// { 1,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.03*c-0.21*s;
    dr = dr+0.38*c-0.06*s;
    db = db-0.02*s;
    u = 2*um3-um5; //{ 2,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.16*c+0.05*s;
    dr = dr-0.18*c-0.31*s;
    db = db+0.01*c;
    u = 2*um3-2*um5;// { 2,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.14*c-2.73*s;
    dr = dr+9.23*c+0.48*s;
    u = 2*um3-3*um5;// { 2,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.07*c-0.55*s;
    dr = dr+1.83*c+0.25*s;
    db = db+0.01*c;
    u = 2*um3-4*um5;// { 2,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.02*c-0.08*s;
    dr = dr+0.25*c+0.06*s;
    u = 3*um3-2*um5; //{ 3,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.01*c-0.07*s;
    dr = dr+0.16*c+0.04*s;
    u = 3*um3-3*um5; //{ 3,-3,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.16*c-0.03*s;
    dr = dr+0.08*c-0.64*s;
    u = 3*um3-4*um5; //{ 3,-4,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.04*c-0.01*s;
    dr = dr+0.03*c-0.17*s;
    //{ perturbations by Saturn }
    u = -um6; //{ 0,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.32*s;
    dr = dr+0.01*c;
    u = um3-um6;// { 1,-1,0 }
    c = cosl(u); s = sinl(u);
    dl = dl-0.08*c-0.41*s;
    dr = dr+0.97*c-0.18*s;
    db = db-0.01*s;
    u = um3-2*um6;// { 1,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.04*c+0.10*s;
    dr = dr-0.23*c+0.10*s;
    u = 2*um3-2*um6;// { 2,-2,0 }
    c = cosl(u); s = sinl(u);
    dl = dl+0.04*c+0.10*s;
    dr = dr-0.35*c+0.13*s;


    //{ difference of Earth-Moon-barycentre and centre of the Earth }
    dl = dl+6.45*sinl(dad)-0.42*sinl(dad-dal)+0.18*sinl(dad+dal)
        +0.17*sinl(dad-um3)-0.06*sinl(dad+um3);
    dr = dr+30.76*cosl(dad)-3.06*cosl(dad-dal)+0.85*cosl(dad+dal)
        -0.58*cosl(dad+um3)+0.57*cosl(dad-um3);
    db = db+0.576*sinl(daf);
    //{ long-periodic perturbations }
    dl = dl+6.40*sinl(PI2*(0.6983+0.0561*ts))
        +1.87*sinl(PI2*(0.5764+0.4174*ts))
        +0.27*sinl(PI2*(0.4189+0.3306*ts))
        +0.20*sinl(PI2*(0.3581+2.4814*ts));
    //{ ecliptic coordinates ([rad],[AU]) }
    *l = PI2*frac(0.7859453+um3/PI2
                 +((6191.2+1.1*ts)*ts+dl)/1296.0e3);
    *r = 1.0001398-0.0000007*ts+dr*1.0e-6;
    *b = SEC_IN_RAD*db;

    return;
}


void get_sun_celestial_position(long double tdb, long double coordinates[3])
{
    long double l, b, r;
    long double eps;
    long double se[3], Rx[3][3], precession[3][3], m_pt[3][3], R[3][3];

    get_sun_ecliptic_position(tdb, &l, &b, &r);
    spherical_to_cartesian(l, b, r, se);

    eps = get_eps_mean(tdb);

    rotate_by_x(-eps, Rx);
    get_precession_matrix(tdb, precession);
    transpose(m_pt, precession, 3, 3);

    mult_matrix(m_pt, Rx, R);

    mult_matrix_by_vector(R, se, coordinates);

    return;
}