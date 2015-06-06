//
// Created by user on 29.05.2015.
//

#include <math.h>
#include "moon.h"
#include "nutation.h"
#include "constants.h"
#include "coordinates_converters.h"
#include "rotation_matrix.h"
#include "matrix_operations.h"
#include "precession.h"


void get_moon_celestial_positionl(long double tdb, long double coordinates[3])
{
    long double l, b, r;
    long double xe[3];
    long double eps;
    long double Rx[3][3];
    long double precession[3][3], m_pt[3][3], R[3][3];

    get_moon_ecliptic_positionl(tdb, &l, &b, &r);
    spherical_to_cartesianl(l, b, r, xe);

    eps = get_eps_meanl(tdb);

    rotxl(-eps, Rx);

    get_precession_matrixl(tdb, precession);
    transposel(m_pt, precession, 3, 3);

    mult_matricesl(m_pt, Rx, R);
    mult_matrix_by_vectorl(R, xe, coordinates);

    return;
}

void get_moon_celestial_position(double tdb, double coordinates[3])
{
    double l, b, r;
    double xe[3];
    double eps;
    double Rx[3][3];
    double precession[3][3], m_pt[3][3], R[3][3];

    get_moon_ecliptic_position(tdb, &l, &b, &r);
    spherical_to_cartesian(l, b, r, xe);

    eps = get_eps_mean(tdb);

    rotx(-eps, Rx);

    get_precession_matrix(tdb, precession);
    transpose(m_pt, precession, 3, 3);

    mult_matrices(m_pt, Rx, R);
    mult_matrix_by_vector(R, xe, coordinates);

    return;
}


void get_moon_ecliptic_positionl(long double tdb, long double *l, long double *b, long double *r)
{
    long double cfa[5]; // скорректированные фундаментальные аргументы
    long double ar, al, ab;
//    long double s, c, u;
//    long double al, ag, az, ap, an;
//    long double fl1, fl2, fl3, fl4;
//    long double fp1, fp2, fp3;
//    long double ff1, ff2, ff3, ff4;
//    long double ufs;

//    long double dt = (tdb - MJD2000)/JULIAN_C;
//
//    long double dgam =  -3332.0e-9*sinl(PI2*(0.59734-5.37261*dt))
//                         -539.0e-9*sinl(PI2*(0.35498-5.37899*dt))
//                         -64.0e-9*sinl(PI2*(0.39943-5.37511*dt));
//
//
//    fl1 = 1.000002208; // { factor for terms with variable l power 1 }
//    fl2 = fl1*fl1;     // { factor for terms with variable l power 2 }
//    fl3 = fl2*fl1;     // { factor for terms with variable l power 3 }
//    fl4 = fl3*fl1;     // { factor for terms with variable l power 4 }
//
//    fp1 = 0.997504612-0.002495388*dt; // { factor for terms with var. l' power 1 }
//    fp2 = fp1*fp1;     // { factor for terms with variable l' power 2 }
//    fp3 = fp2*fp1;     // { factor for terms with variable l' power 3 }
//
//    ff1 = 1.000002708+139.978*dgam; // { factor for terms with variable F power 1 }
//    ff2 = ff1*ff1;     // { factor for terms with variable F power 2 }
//    ff3 = ff2*ff1;     // { factor for terms with variable F power 3 }
//    ff4 = ff3*ff1;     // { factor for terms with variable F power 4 }
//
    get_corr_fund_argsl(tdb, cfa);
//
//    //{ Solar perturbations }
//
//    al = 0.0;            az = 0.0;
//    ag = 0.0;         ap = 3422.700;     ////{ l l'F D }
//        u = 4*cfa[4];
//        s = sinl(u); c = cosl(u);                         ////{ 0 0 0 4 }
//        al = al+13.902*s;    az = az+14.06*s;
//        ag = ag-0.001*c;  ap = ap+0.2607*c;  ////{ 0 0 0 4 }
//        u = 3*cfa[4];
//        s = sinl(u); c = cosl(u);                         ////{ 0 0 0 3 }
//        al = al+0.403*s;     az = az-4.01*s;
//        ag = ag+0.394*c;  ap = ap+0.0023*c;  ////{ 0 0 0 3 }
//        u = 2*cfa[4];
//        s = sinl(u); c = cosl(u);                        // //{ 0 0 0 2 }
//        al = al+2369.912*s;  az = az+2373.36*s;
//        ag = ag+0.601*c;  ap = ap+28.2333*c; ////{ 0 0 0 2 }
//        u = cfa[4];
//        s = sinl(u); c = cosl(u);                        // //{ 0 0 0 1 }
//        al = al-125.154*s;   az = az-112.79*s;
//        ag = ag-0.725*c;  ap = ap-0.9781*c;  ////{ 0 0 0 1 }
//        u = cfa[1]+4*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                         ////{ 1 0 0 4 }
//        al = al+1.979*s;     az = az+6.98*s;
//        ag = ag-0.445*c;  ap = ap+0.0433*c; // //{ 1 0 0 4 }
//        u = cfa[1]+2*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                        // //{ 1 0 0 2 }
//        al = al+191.953*s;   az = az+192.72*s;
//        ag = ag+0.029*c;  ap = ap+3.0861*c; // //{ 1 0 0 2 }
//        u = cfa[1]+cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                      //   //{ 1 0 0 1 }
//        al = al-8.466*s;     az = az-13.51*s;
//        ag = ag+0.455*c;  ap = ap-0.1093*c; // //{ 1 0 0 1 }
//        u = cfa[1];
//        s = fl1*sinl(u); c = fl1*cosl(u);                         ////{ 1 0 0 0 }
//        al = al+22639.500*s; az = az+22609.07*s;
//        ag = ag+0.079*c;  ap = ap+186.5398*c;// //{1 0 0 0 }
//        u = cfa[1]-cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                       //  //{ 1 0 0-1 }
//        al = al+18.609*s;    az = az+3.59*s;
//        ag = ag-0.094*c;  ap = ap+0.0118*c;  ////{ 1 0 0-1 }
//        u = cfa[1]-2*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                        // //{ 1 0 0-2 }
//        al = al-4586.465*s;  az = az-4578.13*s;
//        ag = ag-0.077*c;  ap = ap+34.3117*c;// //{ 1 0 0-2 }
//        u = cfa[1]-3*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                       //  //{ 1 0 0-3 }
//        al = al+3.215*s;     az = az+5.44*s;
//        ag = ag+0.192*c;  ap = ap-0.0386*c; // //{ 1 0 0-3 }
//        u = cfa[1]-4*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                      //   //{ 1 0 0-4 }
//        al = al-38.428*s;    az = az-38.64*s;
//        ag = ag+0.001*c;  ap = ap+0.6008*c; // //{ 1 0 0-4 }
//        u = cfa[1]-6*cfa[4];
//        s = fl1*sinl(u); c = fl1*cosl(u);                        // //{ 1 0 0-6 }
//        al = al-0.393*s;     az = az-1.43*s;
//        ag = ag-0.092*c;  ap = ap+0.0086*c;  ////{ 1 0 0-6 }
//        u = cfa[2]+4*cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al-0.289*s;     az = az-1.59*s;
//        ag = ag+0.123*c;  ap = ap-0.0053*c; // //{ 0 1 0 4 }
//        u = cfa[2]+2*cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al-24.420*s;    az = az-25.10*s;
//        ag = ag+0.040*c;  ap = ap-0.3000*c;  ////{ 0 1 0 2 }
//        u = cfa[2]+cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al+18.023*s;    az = az+17.93*s;
//        ag = ag+0.007*c;  ap = ap+0.1494*c; // //{ 0 1 0 1 }
//        u = cfa[2];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al-668.146*s;   az = az-126.98*s;
//        ag = ag-1.302*c;  ap = ap-0.3997*c; // //{ 0 1 0 0 }
//        u = cfa[2]-cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al+0.560*s;     az = az+0.32*s;
//        ag = ag-0.001*c;  ap = ap-0.0037*c; // //{ 0 1 0-1 }
//        u = cfa[2]-2*cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al-165.145*s;   az = az-165.06*s;
//        ag = ag+0.054*c;  ap = ap+1.9178*c;  ////{ 0 1 0-2 }
//        u = cfa[2]-4*cfa[4];
//        s = fp1*sinl(u); c = fp1*cosl(u);
//        al = al-1.877*s;     az = az-6.46*s;
//        ag = ag-0.416*c;  ap = ap+0.0339*c; // //{ 0 1 0-4 }
//        u = 2*cfa[1]+4*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al+0.213*s;     az = az+1.02*s;
//        ag = ag-0.074*c;  ap = ap+0.0054*c;  ////{ 2 0 0 4 }
//        u = 2*cfa[1]+2*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al+14.387*s;    az = az+14.78*s;
//        ag = ag-0.017*c;  ap = ap+0.2833*c;  ////{ 2 0 0 2 }
//        u = 2*cfa[1]+cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al-0.586*s;     az = az-1.20*s;
//        ag = ag+0.054*c;  ap = ap-0.0100*c; // //{ 2 0 0 1 }
//        u = 2*cfa[1];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al+769.016*s;   az = az+767.96*s;
//        ag = ag+0.107*c;  ap = ap+10.1657*c; ////{ 2 0 0 0 }
//        u = 2*cfa[1]-cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al+1.750*s;     az = az+2.01*s;
//        ag = ag-0.018*c;  ap = ap+0.0155*c;  // //{ 2 0 0-1 }
//        u = 2*cfa[1]-2*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al-211.656*s;   az = az-152.53*s;
//        ag = ag+5.679*c;  ap = ap-0.3039*c; // //{ 2 0 0-2 }
//        u = 2*cfa[1]-3*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al+1.225*s;     az = az+0.91*s;
//        ag = ag-0.030*c;  ap = ap-0.0088*c; // //{ 2 0 0-3 }
//        u = 2*cfa[1]-4*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al-30.773*s;    az = az-34.07*s;
//        ag = ag-0.308*c;  ap = ap+0.3722*c;  ////{ 2 0 0-4 }
//        u = 2*cfa[1]-6*cfa[4];
//        s = fl2*sinl(u); c = fl2*cosl(u);
//        al = al-0.570*s;     az = az-1.40*s;
//        ag = ag-0.074*c;  ap = ap+0.0109*c;  ////{ 2 0 0-6 }
//        u = cfa[1]+cfa[2]+2*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);   // //{ 1 1 0 2 }
//        al = al-2.921*s;     az = az-11.75*s;
//        ag = ag+0.787*c;  ap = ap-0.0484*c;
//        u = cfa[1]+cfa[2]+cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);  //  //{ 1 1 0 1 }
//        al = al+1.267*s;     az = az+1.52*s;
//        ag = ag-0.022*c;  ap = ap+0.0164*c;
//        u = cfa[1]+cfa[2];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al-109.673*s;   az = az-115.18*s;
//        ag = ag+0.461*c;  ap = ap-0.9490*c; // //{ 1 1 0 0 }
//        u = cfa[1]+cfa[2]-2*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al-205.962*s;   az = az-182.36*s;
//        ag = ag+2.056*c;  ap = ap+1.4437*c; // //{ 1 1 0-2 }
//        u = cfa[1]+cfa[2]-3*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+0.233*s;     az = az+0.36*s;
//        ag = ag+0.012*c;  ap = ap-0.0025*c; // //{ 1 1 0-3 }
//        u = cfa[1]+cfa[2]-4*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al-4.391*s;     az = az-9.66*s;
//        ag = ag-0.471*c;  ap = ap+0.0673*c;  ////{ 1 1 0-4 }
//        u = cfa[1]-cfa[2]+4*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+0.283*s;     az = az+1.53*s;
//        ag = ag-0.111*c;  ap = ap+0.0060*c;  ////{ 1-1 0+4 }
//        u = cfa[1]-cfa[2]+2*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+14.577*s;    az = az+31.70*s;
//        ag = ag-1.540*c;  ap = ap+0.2302*c;  //{ 1-1 0 2 }
//        u = cfa[1]-cfa[2];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+147.687*s;   az = az+138.76*s;
//        ag = ag+0.679*c;  ap = ap+1.1528*c;  //{ 1-1 0 0 }
//        u = cfa[1]-cfa[2]-cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al-1.089*s;     az = az+0.55*s;
//        ag = ag+0.021*c; //{0.0   ;}          //{ 1-1 0-1 }
//        u = cfa[1]-cfa[2]-2*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+28.475*s;    az = az+23.59*s;
//        ag = ag-0.443*c;  ap = ap-0.2257*c;  //{ 1-1 0-2 }
//        u = cfa[1]-cfa[2]-3*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al-0.276*s;     az = az-0.38*s;
//        ag = ag-0.006*c;  ap = ap-0.0036*c;  //{ 1-1 0-3 }
//        u = cfa[1]-cfa[2]-4*cfa[4];
//        s = fl1*fp1*sinl(u); c = fl1*fp1*cosl(u);
//        al = al+0.636*s;     az = az+2.27*s;
//        ag = ag+0.146*c;  ap = ap-0.0102*c;  //{ 1-1 0-4 }
//        u = 2*cfa[2]+2*cfa[4];
//        s = fp2*sinl(u); c = fp2*cosl(u);
//        al = al-0.189*s;     az = az-1.68*s;
//        ag = ag+0.131*c;  ap = ap-0.0028*c;  //{ 0 2 0 2 }
//        u = 2*cfa[2];
//        s = fp2*sinl(u); c = fp2*cosl(u);
//        al = al-7.486*s;     az = az-0.66*s;
//        ag = ag-0.037*c;  ap = ap-0.0086*c;  //{ 0 2 0 0 }
//        u = 2*cfa[2]-2*cfa[4];
//        s = fp2*sinl(u); c = fp2*cosl(u);
//        al = al-8.096*s;     az = az-16.35*s;
//        ag = ag-0.740*c;  ap = ap+0.0918*c;  //{ 0 2 0-2 }
//        u = 2*cfa[3]+2*cfa[4];
//        s = ff2*sinl(u); c = ff2*cosl(u);
//        al = al-5.741*s;     az = az-0.04*s;
//        //{ 0.0  ;         }
//        ap = ap-0.0009*c;  //{ 0 0 2 2 }
//        u = 2*cfa[3]+cfa[4];
//        s = ff2*sinl(u);       //{ 0 0 2 1 }
//        al = al+0.255*s;   //{ 0.0 ;
////            0.0  ;           0.0   ; }
//        u = 2*cfa[3];
//        s = ff2*sinl(u); c = ff2*cosl(u);        //{ 0 0 2 0 }
//        al = al-411.608*s;   az = az-0.20*s;
//        //{ 0.0  ;         }
//        ap = ap-0.0124*c;
//        u = 2*cfa[3]-cfa[4];
//        s = ff2*sinl(u); c = ff2*cosl(u);
//        al = al+0.584*s;     az = az+0.84*s;
//        //{ 0.0  ;         }
//        ap = ap+0.0071*c;  //{ 0 0 2-1 }
//        u = 2*cfa[3]-2*cfa[4];
//        s = ff2*sinl(u); c = ff2*cosl(u);
//        al = al-55.173*s;    az = az-52.14*s;
//        //{ 0.0  ;         }
//        ap = ap-0.1052*c;  //{ 0 0 2-2 }
//        u = 2*cfa[3]-3*cfa[4];
//        s = ff2*sinl(u); c = ff2*cosl(u);
//        al = al+0.254*s;     az = az+0.25*s;
//        //{ 0.0  ;         }
//        ap = ap-0.0017*c;  //{ 0 0 2-3 }
//        u = 2*cfa[3]-4*cfa[4];
//        s = ff2*sinl(u); c = ff2*cosl(u);
//        al = al+0.025*s;     az = az-1.67*s;
//        //{ 0.0  ;         }
//        ap = ap+0.0031*c;  //{ 0 0 2-4 }
//        u = 3*cfa[1]+2*cfa[4];
//        s = fl3*sinl(u); c = fl3*cosl(u);
//        al = al+1.060*s;     az = az+2.96*s;
//        ag = ag-0.166*c;  ap = ap+0.0243*c;  //{ 3 0 0+2 }
//        u = 3*cfa[1];
//        s = fl3*sinl(u); c = fl3*cosl(u);
//        al = al+36.124*s;    az = az+50.64*s;
//        ag = ag-1.300*c;  ap = ap+0.6215*c;  //{ 3 0 0 0 }
//        u = 3*cfa[1]-2*cfa[4];
//        s = fl3*sinl(u); c = fl3*cosl(u);
//        al = al-13.193*s;    az = az-16.40*s;
//        ag = ag+0.258*c;  ap = ap-0.1187*c;  //{ 3 0 0-2 }
//        u = 3*cfa[1]-4*cfa[4];
//        s = fl3*sinl(u); c = fl3*cosl(u);
//        al = al-1.187*s;     az = az-0.74*s;
//        ag = ag+0.042*c;  ap = ap+0.0074*c;  //{ 3 0 0-4 }
//        u = 3*cfa[1]-6*cfa[4];
//        s = fl3*sinl(u); c = fl3*cosl(u);
//        al = al-0.293*s;     az = az-0.31*s;
//        ag = ag-0.002*c;  ap = ap+0.0046*c;  //{ 3 0 0-6 }
//        u = 2*cfa[1]+cfa[2]+2*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-0.290*s;     az = az-1.45*s;
//        ag = ag+0.116*c;  ap = ap-0.0051*c;  //{ 2 1 0 2 }
//        u = 2*cfa[1]+cfa[2];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-7.649*s;     az = az-10.56*s;
//        ag = ag+0.259*c;  ap = ap-0.1038*c;  //{ 2 1 0 0 }
//        u = 2*cfa[1]+cfa[2]-2*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-8.627*s;     az = az-7.59*s;
//        ag = ag+0.078*c;  ap = ap-0.0192*c;  //{ 2 1 0-2 }
//        u = 2*cfa[1]+cfa[2]-4*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-2.740*s;     az = az-2.54*s;
//        ag = ag+0.022*c;  ap = ap+0.0324*c;  //{ 2 1 0-4 }
//        u = 2*cfa[1]-cfa[2]+2*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al+1.181*s;     az = az+3.32*s;
//        ag = ag-0.212*c;  ap = ap+0.0213*c;  //{ 2-1 0+2 }
//        u = 2*cfa[1]-cfa[2];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al+9.703*s;     az = az+11.67*s;
//        ag = ag-0.151*c;  ap = ap+0.1268*c;  //{ 2-1 0 0 }
//        u = 2*cfa[1]-cfa[2]-cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-0.352*s;     az = az-0.37*s;
//        ag = ag+0.001*c;  ap = ap-0.0028*c;  //{ 2-1 0-1 }
//        u = 2*cfa[1]-cfa[2]-2*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al-2.494*s;     az = az-1.17*s;
//        ag = ag-0.003*c;  ap = ap-0.0017*c;  //{ 2-1 0-2 }
//        u = 2*cfa[1]-cfa[2]-4*cfa[4];
//        s = fl2*fp1*sinl(u); c = fl2*fp1*cosl(u);
//        al = al+0.360*s;     az = az+0.20*s;
//        ag = ag-0.012*c;  ap = ap-0.0043*c;  //{ 2-1 0-4 }
//        u = cfa[1]+2*cfa[2];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al-1.167*s;     az = az-1.25*s;
//        ag = ag+0.008*c;  ap = ap-0.0106*c;  //{ 1 2 0 0 }
//        u = cfa[1]+2*cfa[2]-2*cfa[4];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al-7.412*s;     az = az-6.12*s;
//        ag = ag+0.117*c;  ap = ap+0.0484*c;  //{ 1 2 0-2 }
//        u = cfa[1]+2*cfa[2]-4*cfa[4];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al-0.311*s;     az = az-0.65*s;
//        ag = ag-0.032*c;  ap = ap+0.0044*c;  //{ 1 2 0-4 }
//        u = cfa[1]-2*cfa[2]+2*cfa[4];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al+0.757*s;     az = az+1.82*s;
//        ag = ag-0.105*c;  ap = ap+0.0112*c;  //{ 1-2 0 2 }
//        u = cfa[1]-2*cfa[2];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al+2.580*s;     az = az+2.32*s;
//        ag = ag+0.027*c;  ap = ap+0.0196*c;  //{ 1-2 0 0 }
//        u = cfa[1]-2*cfa[2]-2*cfa[4];
//        s = fl1*fp2*sinl(u); c = fl1*fp2*cosl(u);
//        al = al+2.533*s;     az = az+2.40*s;
//        ag = ag-0.014*c;  ap = ap-0.0212*c;  //{ 1-2 0-2 }
//        u = 3*cfa[2]-2*cfa[4];
//        s = fp3*sinl(u); c = fp3*cosl(u);
//        al = al-0.344*s;     az = az-0.57*s;
//        ag = ag-0.025*c;  ap = ap+0.0036*c;  //{ 0 3 0-2 }
//        u = cfa[1]+2*cfa[3]+2*cfa[4];
//        s = fl1*ff2*sinl(u);
//        al = al-0.992*s;     az = az-0.02*s;
//        //{ 0.0  ;           0.0   ; }         //{ 1 0 2 2 }
//        u = cfa[1]+2*cfa[3];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al-45.099*s;    az = az-0.02*s;
//        //{ 0.0  ; }
//        ap = ap-0.0010*c;  //{ 1 0 2 0 }
//        u = cfa[1]+2*cfa[3]-2*cfa[4];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al-0.179*s;     az = az-9.52*s;
//        //{ 0.0  ; }
//        ap = ap-0.0833*c;  //{ 1 0 2-2 }
//        u = cfa[1]+2*cfa[3]-4*cfa[4];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al-0.301*s;     az = az-0.33*s;
//        //{ 0.0  ; }
//        ap = ap+0.0014*c;  //{ 1 0 2-4 }
//        u = cfa[1]-2*cfa[3]+2*cfa[4];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al-6.382*s;     az = az-3.37*s;
//        //{ 0.0  ; }
//        ap = ap-0.0481*c;  //{ 1 0-2 2 }
//        u = cfa[1]-2*cfa[3];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al+39.528*s;    az = az+85.13*s;
//        //{ 0.0  ; }
//        ap = ap-0.7136*c;  //{ 1 0-2 0 }
//        u = cfa[1]-2*cfa[3]-2*cfa[4];
//        s = fl1*ff2*sinl(u); c = fl1*ff2*cosl(u);
//        al = al+9.366*s;     az = az+0.71*s;
//        //{ 0.0  ; }
//        ap = ap-0.0112*c;  //{ 1 0-2-2 }
//        u = cfa[1]-2*cfa[3]-4*cfa[4];
//        s = fl1*ff2*sinl(u);
//        al = al+0.202*s;     az = az+0.02*s;
//        //{ 0.0  ;           0.0  }            //{ 1 0-2-4 }
//        u = cfa[2]+2*cfa[3];
//        s = fp1*ff2*sinl(u); c = fp1*ff2*cosl(u);
//        al = al+0.415*s;     az = az+0.10*s;
//        //{ 0.0  ; }
//        ap = ap+0.0013*c;  //{ 0 1 2 0 }
//        u = cfa[2]+2*cfa[3]-2*cfa[4];
//        s = fp1*ff2*sinl(u); c = fp1*ff2*cosl(u);
//        al = al-2.152*s;     az = az-2.26*s;
//        //{ 0.0  ; }
//        ap = ap-0.0066*c;  //{ 0 1 2-2 }
//        u = cfa[2]-2*cfa[3]+2*cfa[4];
//        s = fp1*ff2*sinl(u); c = fp1*ff2*cosl(u);
//        al = al-1.440*s;     az = az-1.30*s;
//        //{ 0.0  ; }
//        ap = ap+0.0014*c;  //{ 0 1-2 2 }
//        u = cfa[2]-2*cfa[3]-2*cfa[4];
//        s = fp1*ff2*sinl(u);
//        al = al+0.384*s;     az = az-0.04*s;
//        //{ 0.0  ;           0.0  }            //{ 0 1-2-2 }
//        u = 4*cfa[1];
//        s = fl4*sinl(u); c = fl4*cosl(u);
//        al = al+1.938*s;     az = az+3.60*s;
//        ag = ag-0.145*c;  ap = ap+0.0401*c;  //{ 4 0 0 0 }
//        u = 4*cfa[1]-2*cfa[4];
//        s = fl4*sinl(u); c = fl4*cosl(u);
//        al = al-0.952*s;     az = az-1.58*s;
//        ag = ag+0.052*c;  ap = ap-0.0130*c;  //{ 4 0 0-2 }
//        u = 3*cfa[1]+cfa[2];
//        s = fl3*fp1*sinl(u); c = fl3*fp1*cosl(u);
//        al = al-0.551*s;     az = az-0.94*s;
//        ag = ag+0.032*c;  ap = ap-0.0097*c;  //{ 3 1 0 0 }
//        u = 3*cfa[1]+cfa[2]-2*cfa[4];
//        s = fl3*fp1*sinl(u); c = fl3*fp1*cosl(u);
//        al = al-0.482*s;     az = az-0.57*s;
//        ag = ag+0.005*c;  ap = ap-0.0045*c;  //{ 3 1 0-2 }
//        u = 3*cfa[1]-cfa[2];
//        s = fl3*fp1*sinl(u); c = fl3*fp1*cosl(u);
//        al = al+0.681*s;     az = az+0.96*s;
//        ag = ag-0.026*c;  ap = ap+0.0115*c;  //{ 3-1 0 0 }
//        u = 2*cfa[1]+2*cfa[2]-2*cfa[4];
//        s = fl2*fp2*sinl(u); c = fl2*fp2*cosl(u);
//        al = al-0.297*s;     az = az-0.27*s;
//        ag = ag+0.002*c;  ap = ap-0.0009*c;  //{ 2 2 0-2 }
//        u = 2*cfa[1]-2*cfa[2]-2*cfa[4];
//        s = fl2*fp2*sinl(u); c = fl2*fp2*cosl(u);
//        al = al+0.254*s;     az = az+0.21*s;
//        ag = ag-0.003*c; //{0.0 }             //{ 2-2 0-2 }
//        u = cfa[1]+3*cfa[2]-2*cfa[4];
//        s = fl1*fp3*sinl(u); c = fl1*fp3*cosl(u);
//        al = al-0.250*s;     az = az-0.22*s;
//        ag = ag+0.004*c;  ap = ap+0.0014*c;  //{ 1 3 0-2 }
//        u = 2*cfa[1]+2*cfa[3];
//        s = fl2*ff2*sinl(u); c = fl2*ff2*cosl(u);
//        al = al-3.996*s;   //{ 0.0 ;0.0  ; }
//        ap = ap+0.0004*c;  //{ 2 0 2 0 }
//        u = 2*cfa[1]+2*cfa[3]-2*cfa[4];
//        s = fl2*ff2*sinl(u); c = fl2*ff2*cosl(u);
//        al = al+0.557*s;     az = az-0.75*s;
//        //{ 0.0  ; }
//        ap = ap-0.0090*c;  //{ 2 0 2-2 }
//        u = 2*cfa[1]-2*cfa[3]+2*cfa[4];
//        s = fl2*ff2*sinl(u); c = fl2*ff2*cosl(u);
//        al = al-0.459*s;     az = az-0.38*s;
//        //{ 0.0  ; }
//        ap = ap-0.0053*c;  //{ 2 0-2 2 }
//        u = 2*cfa[1]-2*cfa[3];
//        s = fl2*ff2*sinl(u); c = fl2*ff2*cosl(u);
//        al = al-1.298*s;     az = az+0.74*s;
//        //{ 0.0  ; }
//        ap = ap+0.0004*c;  //{ 2 0-2 0 }
//        u = 2*cfa[1]-2*cfa[3]-2*cfa[4];
//        s = fl2*ff2*sinl(u); c = fl2*ff2*cosl(u);
//        al = al+0.538*s;     az = az+1.14*s;
//        //{ 0.0  ; }
//        ap = ap-0.0141*c;  //{ 2 0-2-2 }
//        u = cfa[1]+cfa[2]+2*cfa[3];
//        s = fl1*fp1*ff2*sinl(u);
//        al = al+0.263*s;     az = az+0.02*s;
//        //{ 0.0  ;           0.0 }             //{ 1 1 2 0 }
//        u = cfa[1]+cfa[2]-2*cfa[3]-2*cfa[4];
//        s = fl1*fp1*ff2*sinl(u); c = fl1*fp1*ff2*cosl(u);
//        al = al+0.426*s;     az = az+0.07*s;
//        //{ 0.0  ; }
//        ap = ap-0.0006*c;  //{ 1 1-2-2 }
//        u = cfa[1]-cfa[2]+2*cfa[3];
//        s = fl1*fp1*ff2*sinl(u); c = fl1*fp1*ff2*cosl(u);
//        al = al-0.304*s;     az = az+0.03*s;
//        //{ 0.0  ; }
//        ap = ap+0.0003*c;  //{ 1-1 2 0 }
//        u = cfa[1]-cfa[2]-2*cfa[3]+2*cfa[4];
//        s = fl1*fp1*ff2*sinl(u); c = fl1*fp1*ff2*cosl(u);
//        al = al-0.372*s;     az = az-0.19*s;
//        //{ 0.0  ; }
//        ap = ap-0.0027*c;  //{ 1-1-2 2 }
//        u = 4*cfa[3];
//        s = ff4*sinl(u);
//        al = al+0.418*s;   //{ 0.0 ;
//           // 0.0  ;          0.0   ; }            //{ 0 0 4 0 }
//        u = 3*cfa[1]+2*cfa[3];
//        s = fl3*ff2*sinl(u);
//        al = al-0.330*s;     az = az-0.04*s;
//        //{ 0.0  ;   0.0   ;}                   //{ 3 0 2 0 }
//
//        //{ Solar perturbations in latitude }
//        u = cfa[3]-2*cfa[4];     s = ff1*sinl(u);
//        an = -526.069*s;  //{-526.069  0 0 1-2 }
//        u = cfa[3]-4*cfa[4];     s = ff1*sinl(u);
//        an = an-3.352*s;  //{  -3.352  0 0 1-4 }
//        u = cfa[1]+cfa[3]-2*cfa[4]; s = fl1*ff1*sinl(u);
//        an = an+44.297*s; //{ +44.297 +1 0 1-2 }
//        u = cfa[1]+cfa[3]-4*cfa[4]; s = fl1*ff1*sinl(u);
//        an = an-6.000*s;  //{  -6.000 +1 0 1-4 }
//        u = -cfa[1]+cfa[3];      s = fl1*ff1*sinl(u);
//        an = an+20.599*s; //{ +20.599 -1 0 1 0 }
//        u = -cfa[1]+cfa[3]-2*cfa[4];s = fl1*ff1*sinl(u);
//        an = an-30.598*s; //{ -30.598 -1 0 1-2 }
//        u = -2*cfa[1]+cfa[3];    s = fl2*ff1*sinl(u);
//        an = an-24.649*s; //{ -24.649 -2 0 1 0 }
//        u = -2*cfa[1]+cfa[3]-2*cfa[4]; s = fl2*ff1*sinl(u);
//        an = an-2.0*s; //{  -2.000 -2 0 1-2 }
//        u = cfa[2]+cfa[3]-2*cfa[4]; s = fp1*ff1*sinl(u);
//        an = an-22.571*s; //{ -22.571  0+1 1-2 }
//        u = -cfa[2]+cfa[3]-2*cfa[4];s = fp1*ff1*sinl(u);
//        an = an+10.985*s; //{ +10.985  0-1 1-2 }
//
//        //{ planetary: perturbations
//            //of ecliptic latitude by Venus and Jupiter }
//        al = al+0.82*sinl(PI2*(0.7736-62.5512*dt))
//            +0.31*sinl(PI2*(0.0466-125.1025*dt))
//            +0.35*sinl(PI2*(0.5785-25.1042*dt))
//            +0.66*sinl(PI2*(0.4591+1335.8075*dt))
//            +0.64*sinl(PI2*(0.3130-91.5680*dt))
//            +1.14*sinl(PI2*(0.1480+1331.2898*dt))
//            +0.21*sinl(PI2*(0.5918+1056.5859*dt))
//            +0.44*sinl(PI2*(0.5784+1322.8595*dt))
//            +0.24*sinl(PI2*(0.2275-5.7374*dt))
//            +0.28*sinl(PI2*(0.2965+2.6929*dt))
//            +0.33*sinl(PI2*(0.3132+6.3368*dt));
//
//        *l = cfa[0]+SEC_IN_RAD*al; //{ ecliptic longitude of the Moon in radian }
//        ufs = cfa[3]+SEC_IN_RAD*az; //{ angle S=F+deltaS in radian }
//        *b = SEC_IN_RAD*(ff1*(18519.7+ag)*sinl(ufs)-6.24*sinl(3*ufs)+an); //{ latitude }
//        *r = R0 / (SEC_IN_RAD*0.999953253*ap); //{ range in km }
//
    
    

    ar= //{3422.70}
          +0.260968 * cosl(4*cfa[4])
         +28.233869 * cosl(2*cfa[4])
          +0.043566 * cosl(cfa[1] + 4*cfa[4])
           +3.08589 * cosl(cfa[1] + 2*cfa[4])
        +186.539296 * cosl(cfa[1])
         +34.311569 * cosl(cfa[1] - 2*cfa[4])
           +0.60071 * cosl(cfa[1] - 4*cfa[4])
          -0.300334 * cosl(cfa[2] + 2*cfa[4])
          -0.399822 * cosl(cfa[2])
          +1.916735 * cosl(cfa[2] - 2*cfa[4])
          +0.034671 * cosl(cfa[2] - 4*cfa[4])
          -0.977818 * cosl(cfa[4])
          +0.282799 * cosl(2*cfa[1] + 2*cfa[4])
         +10.165933 * cosl(2*cfa[1])
          -0.304041 * cosl(2*cfa[1] - 2*cfa[4])
          +0.372337 * cosl(2*cfa[1] - 4*cfa[4])
          -0.949147 * cosl(cfa[1] + cfa[2])
          +1.443617 * cosl(cfa[1] + cfa[2] - 2*cfa[4])
          +0.067283 * cosl(cfa[1] + cfa[2] - 4*cfa[4])
          +0.229935 * cosl(cfa[1] - cfa[2] + 2*cfa[4])
          +1.152852 * cosl(cfa[1] - cfa[2])
          -0.225821 * cosl(cfa[1] - cfa[2] - 2*cfa[4])
          -0.008639 * cosl(2*cfa[2])
          +0.091646 * cosl(2*cfa[2] - 2*cfa[4])
          -0.012103 * cosl(2*cfa[3])
          -0.105291 * cosl(2*cfa[3] - 2*cfa[4])
          -0.109456 * cosl(cfa[1] + cfa[4])
          +0.011715 * cosl(cfa[1] - cfa[4])
          -0.038258 * cosl(cfa[1] - 3*cfa[4])
          +0.149444 * cosl(cfa[2] + cfa[4])
          -0.225821 * cosl(cfa[1] - cfa[2] - 2*cfa[4])
          -0.008639 * cosl(2*cfa[2])
          +0.091646 * cosl(2*cfa[2] - 2*cfa[4])
          -0.012103 * cosl(2*cfa[3])
          -0.105291 * cosl(2*cfa[3] - 2*cfa[4])
          -0.109456 * cosl(cfa[1] + cfa[4])
          +0.011715 * cosl(cfa[1] - cfa[4])
          -0.038258 * cosl(cfa[1] - 3*cfa[4])
          -0.118714 * cosl(3*cfa[1] - 2*cfa[4])
          -0.047853 * cosl(cfa[1] - 2*cfa[3] + 2*cfa[4])
          -0.708093 * cosl(cfa[1] - 2*cfa[3])
          -0.048117 * cosl(cfa[1] + cfa[2] + 2*cfa[4])
          +0.621546 * cosl(3*cfa[1])
          -0.103337 * cosl(2*cfa[1] + cfa[2])
          -0.083228 * cosl(cfa[1] + 2*cfa[3] - 2*cfa[4]);

    //  for longitude in arcsec
    al =  +13.90 * sinl(4*cfa[4])
        +2369.92 * sinl(2*cfa[4])
            +1.98 * sinl(cfa[1] + 4*cfa[4])
          +191.96 * sinl(cfa[1] + 2*cfa[4])
         -4586.47 * sinl(cfa[1] - 2*cfa[4])
           -38.43 * sinl(cfa[1] - 4*cfa[4])
           -24.42 * sinl(cfa[2] + 2*cfa[4])
          -668.15 * sinl(cfa[2])
          -165.15 * sinl(cfa[2] - 2*cfa[4])
            -1.88 * sinl(cfa[2] - 4*cfa[4])
          -125.15 * sinl(cfa[4])
           +14.38 * sinl(2*cfa[1] + 2*cfa[4])
          +769.02 * sinl(2*cfa[1])
          -211.66 * sinl(2*cfa[1] - 2*cfa[4])
           -30.77 * sinl(2*cfa[1] - 4*cfa[4])
            -2.92 * sinl(cfa[1] + cfa[2] + 2*cfa[4])
          -109.67 * sinl(cfa[1] + cfa[2])
          -205.96 * sinl(cfa[1] + cfa[2] - 2*cfa[4])
            -4.39 * sinl(cfa[1] + cfa[2] - 4*cfa[4])
           +14.57 * sinl(cfa[1] - cfa[2] + 2*cfa[4])
          +147.69 * sinl(cfa[1] - cfa[2])
           +28.47 * sinl(cfa[1] - cfa[2] - 2*cfa[4])
            -7.49 * sinl(2*cfa[2])
            -8.09 * sinl(2*cfa[2] - 2*cfa[4])
            -5.74 * sinl(2*cfa[3] + 2*cfa[4])
          -411.60 * sinl(2*cfa[3])
           -55.17 * sinl(2*cfa[3] - 2*cfa[4])
            -8.46 * sinl(cfa[1] + cfa[4])
           +18.61 * sinl(cfa[1] - cfa[4])
            +3.21 * sinl(cfa[1] - 3*cfa[4])
           +18.02 * sinl(cfa[2] + cfa[4])
            +0.56 * sinl(cfa[2] - cfa[4])
           +36.12 * sinl(3*cfa[1])
           -13.19 * sinl(3*cfa[1] - 2*cfa[4])
            -7.65 * sinl(2*cfa[1] + cfa[2])
            +9.70 * sinl(2*cfa[1] - cfa[2])
            -2.49 * sinl(2*cfa[1] - cfa[2] - 2*cfa[4])
            -0.99 * sinl(cfa[1] + 2*cfa[3] + 2*cfa[4])
           -45.10 * sinl(cfa[1] + 2*cfa[3])
            -6.38 * sinl(cfa[1] - 2*cfa[3] + 2*cfa[4])
           +39.53 * sinl(cfa[1] - 2*cfa[3])
            +1.75 * sinl(2*cfa[1] - cfa[4])
        +22639.50 * sinl(cfa[1])
            -0.57 * sinl(2*cfa[1] - 6*cfa[4])
            +0.64 * sinl(cfa[1] - cfa[2] - 4*cfa[4])
            +1.06 * sinl(3*cfa[1] + 2*cfa[4])
            -1.19 * sinl(3*cfa[1] - 4*cfa[4])
            -8.63 * sinl(2*cfa[1] + cfa[2] - 2*cfa[4])
            -2.74 * sinl(2*cfa[1] + cfa[2] - 4*cfa[4])
            +1.18 * sinl(2*cfa[1] - cfa[2] + 2*cfa[4])
            -1.17 * sinl(cfa[1] + 2*cfa[2])
            -7.41 * sinl(cfa[1] + 2*cfa[2] - 2*cfa[4])
            +0.76 * sinl(cfa[1] - 2*cfa[2] + 2*cfa[4])
            +2.58 * sinl(cfa[1] - 2*cfa[2])
            +2.53 * sinl(cfa[1] - 2*cfa[2] - 2*cfa[4])
            +9.37 * sinl(cfa[1] - 2*cfa[3] - 2*cfa[4])
            -2.15 * sinl(cfa[2] + 2*cfa[3] - 2*cfa[4])
            -1.44 * sinl(cfa[2] - 2*cfa[3] + 2*cfa[4])
            -0.59 * sinl(2*cfa[1] + cfa[4])
            +1.22 * sinl(2*cfa[1] - 3*cfa[4])
            +1.27 * sinl(cfa[1]+cfa[2] + cfa[4])
            -1.09 * sinl(cfa[1]-cfa[2] - cfa[4])
            +0.58 * sinl(2*cfa[3] - cfa[4])
            +1.94 * sinl(4*cfa[1])
            -0.95 * sinl(4*cfa[1] - 2*cfa[4])
            -0.55 * sinl(3*cfa[1] + cfa[2])
            +0.67 * sinl(3*cfa[1] - cfa[2])
            -4.00 * sinl(2*cfa[1] + 2*cfa[3])
            +0.56 * sinl(2*cfa[1] + 2*cfa[3] - 2*cfa[4])
            -1.30 * sinl(2*cfa[1] - 2*cfa[3])
            +0.54 * sinl(2*cfa[1] - 2*cfa[3] - 2*cfa[4]);

    //  for latitude in arcsec
    ab =     117.26 * sinl(cfa[3] + 2*cfa[4])
          +18461.35 * sinl(cfa[3])
            -623.66 * sinl(cfa[3] - 2*cfa[4])
              -3.67 * sinl(cfa[3] - 4*cfa[4])
             +15.12 * sinl(cfa[1] + cfa[3] + 2*cfa[4])
            -166.58 * sinl(cfa[1] + cfa[3] - 2*cfa[4])
              -6.58 * sinl(cfa[1] + cfa[3] - 4*cfa[4])
              +3.00 * sinl(-cfa[1] + cfa[3] + 4*cfa[4])
            +199.49 * sinl(-cfa[1] + cfa[3] + 2*cfa[4])
            -999.69 * sinl(-cfa[1] + cfa[3])
             -33.36 * sinl(-cfa[1] + cfa[3] - 2*cfa[4])
              -6.48 * sinl(cfa[2] + cfa[3])
             -29.65 * sinl(cfa[2] + cfa[3] - 2*cfa[4])
              +7.98 * sinl(-cfa[2] + cfa[3] + 2*cfa[4])
              +4.86 * sinl(-cfa[2] + cfa[3])
              -5.38 * sinl(cfa[3] + cfa[4])
              +4.81 * sinl(cfa[3] - cfa[4])
             -15.57 * sinl(2*cfa[1] + cfa[3] - 2*cfa[4])
             -31.76 * sinl(-2*cfa[1] + cfa[3])
              -5.33 * sinl(cfa[1] + cfa[2] + cfa[3])
              +8.89 * sinl(-cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
              +6.75 * sinl(cfa[1] - cfa[2] + cfa[3])
              -5.65 * sinl(-cfa[1] + cfa[2] + cfa[3])
              -1.02 * sinl(cfa[1] + 3*cfa[3])
              +1.19 * sinl(cfa[3] + 4*cfa[4])
           +1010.16 * sinl(cfa[1] + cfa[3])
              -1.26 * sinl(cfa[2] + cfa[3] + 2*cfa[4])
             +12.12 * sinl(-cfa[2] + cfa[3] - 2*cfa[4])
              -6.29 * sinl(3*cfa[3])
              -2.18 * sinl(3*cfa[3] - 2*cfa[4])
              +1.51 * sinl(2*cfa[1] + cfa[3] + 2*cfa[4])
             +61.91 * sinl(2*cfa[1] + cfa[3])
              +2.41 * sinl(-2*cfa[1] + cfa[3] + 4*cfa[4])
              -1.62 * sinl(-2*cfa[1] + cfa[3] + 2*cfa[4])
              -2.14 * sinl(-2*cfa[1] + cfa[3] - 2*cfa[4])
              -7.45 * sinl(cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
              +5.08 * sinl(-cfa[1] - cfa[2] + cfa[3])
              +1.13 * sinl(cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
              -1.32 * sinl(-cfa[1] + cfa[2] + cfa[3] + 2*cfa[4])
              -1.77 * sinl(-cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
              -1.09 * sinl(2*cfa[2] + cfa[3] - 2*cfa[4])
              -2.79 * sinl(-cfa[1] + 3*cfa[3])
              +3.98 * sinl(3*cfa[1] + cfa[3])
              -1.51 * sinl(3*cfa[1] + cfa[3] - 2*cfa[4])
              -1.58 * sinl(-3*cfa[1] + cfa[3]);

    *r = R0 / (SEC_IN_RAD*0.999953253 * (3422.7 + ar));
    *l = cfa[0] + SEC_IN_RAD*al;
    *b = SEC_IN_RAD * ab;
}

void get_moon_ecliptic_position(double tdb,double *l, double *b, double *r)
{
    double cfa[5]; // скорректированные фундаментальные аргументы
    double ar, al, ab;

    get_corr_fund_args(tdb, cfa);

    ar= //{3422.70}
              +0.260968 * cos(4*cfa[4])
             +28.233869 * cos(2*cfa[4])
              +0.043566 * cos(cfa[1] + 4*cfa[4])
               +3.08589 * cos(cfa[1] + 2*cfa[4])
            +186.539296 * cos(cfa[1])
             +34.311569 * cos(cfa[1] - 2*cfa[4])
               +0.60071 * cos(cfa[1] - 4*cfa[4])
              -0.300334 * cos(cfa[2] + 2*cfa[4])
              -0.399822 * cos(cfa[2])
              +1.916735 * cos(cfa[2] - 2*cfa[4])
              +0.034671 * cos(cfa[2] - 4*cfa[4])
              -0.977818 * cos(cfa[4])
              +0.282799 * cos(2*cfa[1] + 2*cfa[4])
             +10.165933 * cos(2*cfa[1])
              -0.304041 * cos(2*cfa[1] - 2*cfa[4])
              +0.372337 * cos(2*cfa[1] - 4*cfa[4])
              -0.949147 * cos(cfa[1] + cfa[2])
              +1.443617 * cos(cfa[1] + cfa[2] - 2*cfa[4])
              +0.067283 * cos(cfa[1] + cfa[2] - 4*cfa[4])
              +0.229935 * cos(cfa[1] - cfa[2] + 2*cfa[4])
              +1.152852 * cos(cfa[1] - cfa[2])
              -0.225821 * cos(cfa[1] - cfa[2] - 2*cfa[4])
              -0.008639 * cos(2*cfa[2])
              +0.091646 * cos(2*cfa[2] - 2*cfa[4])
              -0.012103 * cos(2*cfa[3])
              -0.105291 * cos(2*cfa[3] - 2*cfa[4])
              -0.109456 * cos(cfa[1] + cfa[4])
              +0.011715 * cos(cfa[1] - cfa[4])
              -0.038258 * cos(cfa[1] - 3*cfa[4])
              +0.149444 * cos(cfa[2] + cfa[4])
              -0.225821 * cos(cfa[1] - cfa[2] - 2*cfa[4])
              -0.008639 * cos(2*cfa[2])
              +0.091646 * cos(2*cfa[2] - 2*cfa[4])
              -0.012103 * cos(2*cfa[3])
              -0.105291 * cos(2*cfa[3] - 2*cfa[4])
              -0.109456 * cos(cfa[1] + cfa[4])
              +0.011715 * cos(cfa[1] - cfa[4])
              -0.038258 * cos(cfa[1] - 3*cfa[4])
              -0.118714 * cos(3*cfa[1] - 2*cfa[4])
              -0.047853 * cos(cfa[1] - 2*cfa[3] + 2*cfa[4])
              -0.708093 * cos(cfa[1] - 2*cfa[3])
              -0.048117 * cos(cfa[1] + cfa[2] + 2*cfa[4])
              +0.621546 * cos(3*cfa[1])
              -0.103337 * cos(2*cfa[1] + cfa[2])
              -0.083228 * cos(cfa[1] + 2*cfa[3] - 2*cfa[4]);

    //  for longitude in arcsec
    al =  +13.90 * sin(4*cfa[4])
          +2369.92 * sin(2*cfa[4])
          +1.98 * sin(cfa[1] + 4*cfa[4])
          +191.96 * sin(cfa[1] + 2*cfa[4])
          -4586.47 * sin(cfa[1] - 2*cfa[4])
          -38.43 * sin(cfa[1] - 4*cfa[4])
          -24.42 * sin(cfa[2] + 2*cfa[4])
          -668.15 * sin(cfa[2])
          -165.15 * sin(cfa[2] - 2*cfa[4])
          -1.88 * sin(cfa[2] - 4*cfa[4])
          -125.15 * sin(cfa[4])
          +14.38 * sin(2*cfa[1] + 2*cfa[4])
          +769.02 * sin(2*cfa[1])
          -211.66 * sin(2*cfa[1] - 2*cfa[4])
          -30.77 * sin(2*cfa[1] - 4*cfa[4])
          -2.92 * sin(cfa[1] + cfa[2] + 2*cfa[4])
          -109.67 * sin(cfa[1] + cfa[2])
          -205.96 * sin(cfa[1] + cfa[2] - 2*cfa[4])
          -4.39 * sin(cfa[1] + cfa[2] - 4*cfa[4])
          +14.57 * sin(cfa[1] - cfa[2] + 2*cfa[4])
          +147.69 * sin(cfa[1] - cfa[2])
          +28.47 * sin(cfa[1] - cfa[2] - 2*cfa[4])
          -7.49 * sin(2*cfa[2])
          -8.09 * sin(2*cfa[2] - 2*cfa[4])
          -5.74 * sin(2*cfa[3] + 2*cfa[4])
          -411.60 * sin(2*cfa[3])
          -55.17 * sin(2*cfa[3] - 2*cfa[4])
          -8.46 * sin(cfa[1] + cfa[4])
          +18.61 * sin(cfa[1] - cfa[4])
          +3.21 * sin(cfa[1] - 3*cfa[4])
          +18.02 * sin(cfa[2] + cfa[4])
          +0.56 * sin(cfa[2] - cfa[4])
          +36.12 * sin(3*cfa[1])
          -13.19 * sin(3*cfa[1] - 2*cfa[4])
          -7.65 * sin(2*cfa[1] + cfa[2])
          +9.70 * sin(2*cfa[1] - cfa[2])
          -2.49 * sin(2*cfa[1] - cfa[2] - 2*cfa[4])
          -0.99 * sin(cfa[1] + 2*cfa[3] + 2*cfa[4])
          -45.10 * sin(cfa[1] + 2*cfa[3])
          -6.38 * sin(cfa[1] - 2*cfa[3] + 2*cfa[4])
          +39.53 * sin(cfa[1] - 2*cfa[3])
          +1.75 * sin(2*cfa[1] - cfa[4])
          +22639.50 * sin(cfa[1])
          -0.57 * sin(2*cfa[1] - 6*cfa[4])
          +0.64 * sin(cfa[1] - cfa[2] - 4*cfa[4])
          +1.06 * sin(3*cfa[1] + 2*cfa[4])
          -1.19 * sin(3*cfa[1] - 4*cfa[4])
          -8.63 * sin(2*cfa[1] + cfa[2] - 2*cfa[4])
          -2.74 * sin(2*cfa[1] + cfa[2] - 4*cfa[4])
          +1.18 * sin(2*cfa[1] - cfa[2] + 2*cfa[4])
          -1.17 * sin(cfa[1] + 2*cfa[2])
          -7.41 * sin(cfa[1] + 2*cfa[2] - 2*cfa[4])
          +0.76 * sin(cfa[1] - 2*cfa[2] + 2*cfa[4])
          +2.58 * sin(cfa[1] - 2*cfa[2])
          +2.53 * sin(cfa[1] - 2*cfa[2] - 2*cfa[4])
          +9.37 * sin(cfa[1] - 2*cfa[3] - 2*cfa[4])
          -2.15 * sin(cfa[2] + 2*cfa[3] - 2*cfa[4])
          -1.44 * sin(cfa[2] - 2*cfa[3] + 2*cfa[4])
          -0.59 * sin(2*cfa[1] + cfa[4])
          +1.22 * sin(2*cfa[1] - 3*cfa[4])
          +1.27 * sin(cfa[1]+cfa[2] + cfa[4])
          -1.09 * sin(cfa[1]-cfa[2] - cfa[4])
          +0.58 * sin(2*cfa[3] - cfa[4])
          +1.94 * sin(4*cfa[1])
          -0.95 * sin(4*cfa[1] - 2*cfa[4])
          -0.55 * sin(3*cfa[1] + cfa[2])
          +0.67 * sin(3*cfa[1] - cfa[2])
          -4.00 * sin(2*cfa[1] + 2*cfa[3])
          +0.56 * sin(2*cfa[1] + 2*cfa[3] - 2*cfa[4])
          -1.30 * sin(2*cfa[1] - 2*cfa[3])
          +0.54 * sin(2*cfa[1] - 2*cfa[3] - 2*cfa[4]);

    //  for latitude in arcsec
    ab =     117.26 * sin(cfa[3] + 2*cfa[4])
             +18461.35 * sin(cfa[3])
             -623.66 * sin(cfa[3] - 2*cfa[4])
             -3.67 * sin(cfa[3] - 4*cfa[4])
             +15.12 * sin(cfa[1] + cfa[3] + 2*cfa[4])
             -166.58 * sin(cfa[1] + cfa[3] - 2*cfa[4])
             -6.58 * sin(cfa[1] + cfa[3] - 4*cfa[4])
             +3.00 * sin(-cfa[1] + cfa[3] + 4*cfa[4])
             +199.49 * sin(-cfa[1] + cfa[3] + 2*cfa[4])
             -999.69 * sin(-cfa[1] + cfa[3])
             -33.36 * sin(-cfa[1] + cfa[3] - 2*cfa[4])
             -6.48 * sin(cfa[2] + cfa[3])
             -29.65 * sin(cfa[2] + cfa[3] - 2*cfa[4])
             +7.98 * sin(-cfa[2] + cfa[3] + 2*cfa[4])
             +4.86 * sin(-cfa[2] + cfa[3])
             -5.38 * sin(cfa[3] + cfa[4])
             +4.81 * sin(cfa[3] - cfa[4])
             -15.57 * sin(2*cfa[1] + cfa[3] - 2*cfa[4])
             -31.76 * sin(-2*cfa[1] + cfa[3])
             -5.33 * sin(cfa[1] + cfa[2] + cfa[3])
             +8.89 * sin(-cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
             +6.75 * sin(cfa[1] - cfa[2] + cfa[3])
             -5.65 * sin(-cfa[1] + cfa[2] + cfa[3])
             -1.02 * sin(cfa[1] + 3*cfa[3])
             +1.19 * sin(cfa[3] + 4*cfa[4])
             +1010.16 * sin(cfa[1] + cfa[3])
             -1.26 * sin(cfa[2] + cfa[3] + 2*cfa[4])
             +12.12 * sin(-cfa[2] + cfa[3] - 2*cfa[4])
             -6.29 * sin(3*cfa[3])
             -2.18 * sin(3*cfa[3] - 2*cfa[4])
             +1.51 * sin(2*cfa[1] + cfa[3] + 2*cfa[4])
             +61.91 * sin(2*cfa[1] + cfa[3])
             +2.41 * sin(-2*cfa[1] + cfa[3] + 4*cfa[4])
             -1.62 * sin(-2*cfa[1] + cfa[3] + 2*cfa[4])
             -2.14 * sin(-2*cfa[1] + cfa[3] - 2*cfa[4])
             -7.45 * sin(cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
             +5.08 * sin(-cfa[1] - cfa[2] + cfa[3])
             +1.13 * sin(cfa[1] - cfa[2] + cfa[3] + 2*cfa[4])
             -1.32 * sin(-cfa[1] + cfa[2] + cfa[3] + 2*cfa[4])
             -1.77 * sin(-cfa[1] + cfa[2] + cfa[3] - 2*cfa[4])
             -1.09 * sin(2*cfa[2] + cfa[3] - 2*cfa[4])
             -2.79 * sin(-cfa[1] + 3*cfa[3])
             +3.98 * sin(3*cfa[1] + cfa[3])
             -1.51 * sin(3*cfa[1] + cfa[3] - 2*cfa[4])
             -1.58 * sin(-3*cfa[1] + cfa[3]);

    *r = R0 / (SEC_IN_RAD*0.999953253 * (3422.7 + ar));
    *l = cfa[0] + SEC_IN_RAD*al;
    *b = SEC_IN_RAD * ab;
}
