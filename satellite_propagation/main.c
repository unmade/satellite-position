#include <stdio.h>
#include <math.h>

#include "rotation_matrix.h"
#include "date_converters/date_converters.h"
#include "precession_matrix.h"
#include "constants.h"
//#include "matrix_operations.h"

int main() {
    double a = 0.628318530717959;
    double Rx[3][3];

    rotate_by_x(a, Rx);

    for (int i=0; i<3; i++)
    {
        for (int j=0;j<3;j++)
        {
            printf("%2.12f   ", Rx[i][j]);
        }
        printf("\n");
    }

    printf("\n");

    double jd = utc_to_mjd(2036, 2, 29, 5, 45, 0.0);

    printf("jd=%2.12f\n\n", jd);

    int year, month, day, hour, minute;
    double seconds;
    mjd_to_utc(jd, &year, &month, &day, &hour, &minute, &seconds);

    printf("year=%i, month=%i, day=%i, hour=%i, minute=%i, seconds=%f\n\n",
           year, month, day, hour, minute, seconds);


    printf("∆T = %2.10f\n\n", get_deltaT(2006, 7));


    double tt = utc_to_tt(2003, 11, 15, 15, 35, 0.0);
    printf("TT = %2.12f\n\n", tt);


    double tdb = tt_to_tdb(52958.65004843);
    printf("TDB = %2.12f\n\n", tdb);


    double zeta, theta, z;
    calc_precession_parameters(52905.30143730, &zeta, &theta, &z);
    printf("ζ(t) = %2.12f; θ(t) = %2.12f; z(t) = %2.12f\n\n", zeta*RAD_IN_SEC, theta*RAD_IN_SEC, z*RAD_IN_SEC);


    double precession[3][3];
    calc_precession_matrix(zeta, theta, z, precession);
//    printf("Precession matrix is:\n");
//    print_matrix(sizeof(precession), sizeof(precession[0]), &precession);

    return 0;
}