//
// Created by user on 03.06.2015.
//

#include <assert.h>
#include <math.h>
#include "nutation_test.h"
#include "../nutation.h"

void get_eps_meanl_test(void)
{
    long double eps = get_eps_meanl(59000.5L);

    assert(fabsl(fabsl(eps) - 0.40904647276235173L) < 1e-17);

    return;
}

void get_eps_mean_test(void)
{
    double eps = get_eps_mean(59000.5);

    assert(fabs(fabs(eps) - 0.40904647276235173) < 1e-15);

    return;
}


void get_fund_argsl_test(void)
{
    long double fund_args[5];
    get_fund_argsl(59000.5L, fund_args);

    assert(fabsl(fabsl(fund_args[0]) - 3.1682698134603013L) < 1e-15);
    assert(fabsl(fabsl(fund_args[1]) - 6.0659050554749712L) < 1e-15);
    assert(fabsl(fabsl(fund_args[2]) - 2.5510330536033869L) < 1e-15);
    assert(fabsl(fabsl(fund_args[3]) - 1.5936104120010106L) < 1e-15);
    assert(fabsl(fabsl(fund_args[4]) - 1.9560899054571104L) < 1e-15);

    return;
}

void get_fund_args_test(void)
{
    double fund_args[5];
    get_fund_args(59000.5, fund_args);

    assert(fabs(fabs(fund_args[0]) - 3.168269813460137) < 1e-15);
    assert(fabs(fabs(fund_args[1]) - 6.0659050554748655) < 1e-15);
    assert(fabs(fabs(fund_args[2]) - 2.5510330536033843) < 1e-15);
    assert(fabs(fabs(fund_args[3]) - 1.5936104120008014) < 1e-15);
    assert(fabs(fabs(fund_args[4])- 1.9560899054570753) < 1e-15);

    return;
}


void get_corr_fund_argsl_test(void)
{
    long double fund_args[5];
    get_corr_fund_argsl(59000.5L, fund_args);

    assert(fabsl(fabsl(fund_args[0]) - 3.16834959262001L) < 1e-14);
    assert(fabsl(fabsl(fund_args[1]) - 6.0660073559340422L) < 1e-15);
    assert(fabsl(fabsl(fund_args[2]) - 2.5509953168599497L) < 1e-15);
    assert(fabsl(fabsl(fund_args[3]) - 1.5932158110006065L) < 1e-15);
    assert(fabsl(fabsl(fund_args[4]) - 1.956207421360256L) < 1e-14);

    return;
}

void get_corr_fund_args_test(void)
{
    double fund_args[5];
    get_corr_fund_args(59000.5, fund_args);

    assert(fabs(fabs(fund_args[0]) - 3.1683495926198453) < 1e-15);
    assert(fabs(fabs(fund_args[1]) - 6.0660073559339365) < 1e-15);
    assert(fabs(fabs(fund_args[2]) - 2.5509953168599471) < 1e-15);
    assert(fabs(fabs(fund_args[3]) - 1.5932158110003973) < 1e-15);
    assert(fabs(fabs(fund_args[4]) - 1.9562074213602212) < 1e-15);

    return;
}


void get_nutation_parametersl_test(void)
{
    long double d_psi, d_eps;
    get_nutation_parametersl(55000.5L, &d_psi, &d_eps);

    assert(fabsl(fabsl(d_psi) - 6.7717028675264897e-005L) < 1e-18);
    assert(fabsl(fabsl(d_eps) - 2.1366313021580097e-005L) < 1e-18);

    return;
}

void get_nutation_parameters_test(void)
{
    double d_psi, d_eps;
    get_nutation_parameters(55000.5, &d_psi, &d_eps);

    assert(fabs(fabs(d_psi) - 6.7717028675264897e-005) < 1e-15);
    assert(fabs(fabs(d_eps) - 2.1366313021580097e-005) < 1e-15);

    return;
}


void get_nutation_matrixl_test(void)
{
    long double nutation[3][3];
    get_nutation_matrixl(55000.5L, nutation);

    assert(fabsl(fabsl(nutation[0][0]) - 0.99999999770720205L) < 1e-16);
    assert(fabsl(fabsl(nutation[0][1]) - 6.2129737515251128e-005L) < 1e-17);
    assert(fabsl(fabsl(nutation[0][2]) - 2.6934952791557186e-005L) < 1e-17);

    assert(fabsl(fabsl(nutation[1][0]) - 6.2129162000436897e-005L) < 1e-17);
    assert(fabsl(fabsl(nutation[1][1]) - 0.99999999784170612L) < 1e-16);
    assert(fabsl(fabsl(nutation[1][2]) - 2.1367149742978082e-005L) < 1e-17);

    assert(fabsl(fabsl(nutation[2][0]) - 2.693628026882862e-005L) < 1e-17);
    assert(fabsl(fabsl(nutation[2][1]) - 2.1365476247942066e-005L) < 1e-17);
    assert(fabsl(fabsl(nutation[2][2]) - 0.99999999940897666L) < 1e-16);

    return;
}


void get_nutation_matrix_test(void)
{
    double nutation[3][3];
    get_nutation_matrix(55000.5, nutation);

    assert(fabsl(fabsl(nutation[0][0]) - 0.99999999770720205L) < 1e-15);
    assert(fabsl(fabsl(nutation[0][1]) - 6.2129737515251128e-005L) < 1e-15);
    assert(fabsl(fabsl(nutation[0][2]) - 2.6934952791557186e-005L) < 1e-15);

    assert(fabsl(fabsl(nutation[1][0]) - 6.2129162000436897e-005L) < 1e-15);
    assert(fabsl(fabsl(nutation[1][1]) - 0.99999999784170612L) < 1e-15);
    assert(fabsl(fabsl(nutation[1][2]) - 2.1367149742978082e-005L) < 1e-15);

    assert(fabsl(fabsl(nutation[2][0]) - 2.693628026882862e-005L) < 1e-15);
    assert(fabsl(fabsl(nutation[2][1]) - 2.1365476247942066e-005L) < 1e-15);
    assert(fabsl(fabsl(nutation[2][2]) - 0.99999999940897666L) < 1e-15);

    return;
}


void test_nutation(void)
{
    get_eps_meanl_test();
    get_eps_mean_test();

    get_fund_argsl_test();
    get_fund_args_test();
    get_corr_fund_argsl_test();
    get_corr_fund_args_test();

    get_nutation_parametersl_test();
    get_nutation_parameters_test();

    get_nutation_matrixl_test();
    get_nutation_matrix_test();

    return;
}