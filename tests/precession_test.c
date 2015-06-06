//
// Created by user on 03.06.2015.
//

#include <assert.h>
#include "precession_test.h"
#include "constants.h"
#include "precession.h"


void get_precession_parametersl_test(void)
{
    long double zeta, theta, z;
    get_precession_parametersl(52905.30143730L, &zeta, &theta, &z);

    assert(fabsl(fabsl(zeta) - 0.00041656415882215662L) < 1e-19);
    assert(fabsl(fabsl(theta) - 0.00036202705961655421L) < 1e-19);
    assert(fabsl(fabsl(z) - 0.00041656949403514362L) < 1e-19);
}


void get_precession_parameters_test(void)
{
    double zeta, theta, z;
    get_precession_parameters(52905.30143730, &zeta, &theta, &z);

    assert(fabs(fabs(zeta) - 0.0004165641588221557) < 1e-14);
    assert(fabs(fabs(theta) - 0.0003620270596165534) < 1e-14);
    assert(fabs(fabs(z) - 0.0004165694940351427) < 1e-14);
}


void get_precession_matrixl_test(void)
{
    long double precession[3][3];
    get_precession_matrixl(52905.30143730, precession);

    assert(fabsl(fabsl(precession[0][0]) - 0.999999587412394L) < 1e-14);
    assert(fabsl(fabsl(precession[0][1]) - 0.000833133529175L) < 1e-14);
    assert(fabsl(fabsl(precession[0][2]) - 0.000362027020296L) < 1e-14);

    assert(fabsl(fabsl(precession[1][0]) - 0.000833133529175L) < 1e-14);
    assert(fabsl(fabsl(precession[1][1]) - 0.999999652944190L) < 1e-14);
    assert(fabsl(fabsl(precession[1][2]) - 0.000000150809421L) < 1e-14);

    assert(fabsl(fabsl(precession[2][0]) - 0.000362027020297L) < 1e-14);
    assert(fabsl(fabsl(precession[2][1]) - 0.000000150807490L) < 1e-14);
    assert(fabsl(fabsl(precession[2][2]) - 0.999999934468205L) < 1e-14);

    return;
}


void get_precession_matrix_test(void)
{
    double precession[3][3];
    get_precession_matrix(52905.30143730, precession);

    assert(fabsl(fabsl(precession[0][0]) - 0.999999587412394) < 1e-14);
    assert(fabsl(fabsl(precession[0][1]) - 0.000833133529175) < 1e-14);
    assert(fabsl(fabsl(precession[0][2]) - 0.000362027020296) < 1e-14);

    assert(fabsl(fabsl(precession[1][0]) - 0.000833133529175) < 1e-14);
    assert(fabsl(fabsl(precession[1][1]) - 0.999999652944190) < 1e-14);
    assert(fabsl(fabsl(precession[1][2]) - 0.000000150809421) < 1e-14);

    assert(fabsl(fabsl(precession[2][0]) - 0.000362027020297) < 1e-14);
    assert(fabsl(fabsl(precession[2][1]) - 0.000000150807490) < 1e-14);
    assert(fabsl(fabsl(precession[2][2]) - 0.999999934468205) < 1e-14);

    return;
}


void test_precession(void)
{
    get_precession_parametersl_test();
    get_precession_parameters_test();

    get_precession_matrixl_test();
    get_precession_matrix_test();
}