//
// Created by user on 04.06.2015.
//

#include <assert.h>
#include <math.h>
#include "coordinates_converters_test.h"
#include "coordinates_converters.h"
#include "precession.h"
#include "nutation.h"

void get_mean_equ_to_fixed_matrixl_test(void)
{
    long double precession[3][3], m_p[3][3];
    get_precession_matrixl(52905.30143730L, precession);
    get_mean_equ_to_fixed_matrixl(precession, m_p);

    assert(fabsl(fabsl(m_p[0][0]) - 0.99999958741239448L) < 1e-16);
    assert(fabsl(fabsl(m_p[0][1]) - 0.00083313352917745956L) < 1e-18);
    assert(fabsl(fabsl(m_p[0][2]) - 0.00036202702029796172L) < 1e-18);

    assert(fabsl(fabsl(m_p[1][0]) - 0.00083313352917780922L) < 1e-18);
    assert(fabsl(fabsl(m_p[1][1]) - 0.99999965294418969L) < 1e-16);
    assert(fabsl(fabsl(m_p[1][2]) - 1.5080748990430403e-007L) < 1e-18);

    assert(fabsl(fabsl(m_p[2][0]) - 0.00036202702029715714L) < 1e-18);
    assert(fabsl(fabsl(m_p[2][1]) - 1.5080942139556436e-007L) < 1e-18);
    assert(fabsl(fabsl(m_p[2][2]) - 0.99999993446820479L) < 1e-16);

    return;
}

void get_mean_equ_to_fixed_matrix_test(void)
{
    double precession[3][3], m_p[3][3];
    get_precession_matrix(52905.30143730, precession);
    get_mean_equ_to_fixed_matrix(precession, m_p);

    assert(fabs(fabs(m_p[0][0]) - 0.99999958741239448) < 1e-15);
    assert(fabs(fabs(m_p[0][1]) - 0.00083313352917745772) < 1e-15);
    assert(fabs(fabs(m_p[0][2]) - 0.00036202702029796097) < 1e-15);

    assert(fabs(fabs(m_p[1][0]) - 0.00083313352917780738) < 1e-15);
    assert(fabs(fabs(m_p[1][1]) - 0.99999965294418969) < 1e-15);
    assert(fabs(fabs(m_p[1][2]) - 1.5080748990430337e-007) < 1e-15);

    assert(fabs(fabs(m_p[2][0]) - 0.00036202702029715633) < 1e-15);
    assert(fabs(fabs(m_p[2][1]) - 1.508094213955637e-007) < 1e-15);
    assert(fabs(fabs(m_p[2][2]) - 0.99999993446820479) < 1e-15);

    return;
}


void get_fixed_to_true_equ_matrixl_test(void)
{
    long double precession[3][3], nutation[3][3], m_np[3][3];
    long double tdb = 52905.30143730L;
    get_precession_matrixl(tdb, precession);
    get_nutation_matrixl(tdb, nutation);
    get_fixed_to_true_equ_matrixl(precession, nutation, m_np);

    assert(fabsl(fabsl(m_np[0][0]) - 0.99999964478413361L) < 1e-16);
    assert(fabsl(fabsl(m_np[0][1]) - 0.00077301948693384449L) < 1e-18);
    assert(fabsl(fabsl(m_np[0][2]) - 0.0003359649971987459L) < 1e-18);

    assert(fabsl(fabsl(m_np[1][0]) - 0.00077300948095617849L) < 1e-18);
    assert(fabsl(fabsl(m_np[1][1]) - 0.9999997007807736L) < 1e-16);
    assert(fabsl(fabsl(m_np[1][2]) - 2.9911629751310692e-005L) < 1e-18);

    assert(fabsl(fabsl(m_np[2][0]) - 0.00033598801894424305L) < 1e-18);
    assert(fabsl(fabsl(m_np[2][1]) - 2.9651914998121167e-005L) < 1e-18);
    assert(fabsl(fabsl(m_np[2][2]) - 0.99999994311640594L) < 1e-16);
    
    return;
}

void get_fixed_to_true_equ_matrix_test(void)
{
    double precession[3][3], nutation[3][3], m_np[3][3];
    double tdb = 52905.30143730;
    get_precession_matrix(tdb, precession);
    get_nutation_matrix(tdb, nutation);
    get_fixed_to_true_equ_matrix(precession, nutation, m_np);

    assert(fabs(fabs(m_np[0][0]) - 0.99999964478413372) < 1e-15);
    assert(fabs(fabs(m_np[0][1]) - 0.0007730194869338433) < 1e-15);
    assert(fabs(fabs(m_np[0][2]) - 0.00033596499719874536) < 1e-15);

    assert(fabs(fabs(m_np[1][0]) - 0.00077300948095617719) < 1e-15);
    assert(fabs(fabs(m_np[1][1]) - 0.99999970078077349) < 1e-15);
    assert(fabs(fabs(m_np[1][2]) - 2.9911629751322076e-005) < 1e-15);

    assert(fabs(fabs(m_np[2][0]) - 0.00033598801894424256) < 1e-15);
    assert(fabs(fabs(m_np[2][1]) - 2.9651914998085845e-005) < 1e-15);
    assert(fabs(fabs(m_np[2][2]) - 0.99999994311640583) < 1e-15);
    
    return;
}


void get_true_equ_to_fixed_matrixl_test(void)
{
    long double precession[3][3], nutation[3][3], m_np[3][3], m_pn[3][3];
    long double tdb = 52905.30143730L;
    get_precession_matrixl(tdb, precession);
    get_nutation_matrixl(tdb, nutation);
    get_fixed_to_true_equ_matrixl(precession, nutation, m_np);
    get_true_equ_to_fixed_matrixl(m_np, m_pn);

    assert(fabsl(fabsl(m_pn[0][0]) - 0.99999964478413361L) < 1e-16);
    assert(fabsl(fabsl(m_pn[0][1]) - 0.00077300948095617849L) < 1e-18);
    assert(fabsl(fabsl(m_pn[0][2]) - 0.00033598801894424305L) < 1e-18);
    
    assert(fabsl(fabsl(m_pn[1][0]) - 0.00077301948693384449L) < 1e-18);
    assert(fabsl(fabsl(m_pn[1][1]) - 0.9999997007807736L) < 1e-16);
    assert(fabsl(fabsl(m_pn[1][2]) - 2.9651914998121167e-005L) < 1e-18);
    
    assert(fabsl(fabsl(m_pn[2][0]) - 0.0003359649971987459L) < 1e-18);
    assert(fabsl(fabsl(m_pn[2][1]) - 2.9911629751310692e-005L) < 1e-18);
    assert(fabsl(fabsl(m_pn[2][2]) - 0.99999994311640594L) < 1e-16);

    return;
}

void get_true_equ_to_fixed_matrix_test(void)
{
    double precession[3][3], nutation[3][3], m_np[3][3], m_pn[3][3];
    double tdb = 52905.30143730;
    get_precession_matrix(tdb, precession);
    get_nutation_matrix(tdb, nutation);
    get_fixed_to_true_equ_matrix(precession, nutation, m_np);
    get_true_equ_to_fixed_matrix(m_np, m_pn);

    assert(fabs(fabs(m_pn[0][0]) - 0.99999964478413372) < 1e-15);
    assert(fabs(fabs(m_pn[0][1]) - 0.00077300948095617719) < 1e-15);
    assert(fabs(fabs(m_pn[0][2]) - 0.00033598801894424256) < 1e-15);

    assert(fabs(fabs(m_pn[1][0]) - 0.0007730194869338433) < 1e-15);
    assert(fabs(fabs(m_pn[1][1]) - 0.99999970078077349) < 1e-15);
    assert(fabs(fabs(m_pn[1][2]) - 2.9651914998085845e-005) < 1e-15);

    assert(fabs(fabs(m_pn[2][0]) - 0.00033596499719874536) < 1e-15);
    assert(fabs(fabs(m_pn[2][1]) - 2.9911629751322076e-005) < 1e-15);
    assert(fabs(fabs(m_pn[2][2]) - 0.99999994311640583) < 1e-15);

    return;
}


void get_celes_to_terra_matrixl_test(void)
{
    long double m_ct[3][3];
    get_celes_to_terra_matrixl(51519.114583333336L, 0.38415L, m_ct);

    assert(fabsl(fabsl(m_ct[0][0]) - 0.44919495454806802L) < 1e-16);
    assert(fabsl(fabsl(m_ct[0][1]) - 0.89343376514989092L) < 1e-16);
    assert(fabsl(fabsl(m_ct[0][2]) - 9.932203439718028e-006L) < 1e-16);

    assert(fabsl(fabsl(m_ct[1][0]) - 0.8934337644973771L) < 1e-16);
    assert(fabsl(fabsl(m_ct[1][1]) - 0.44919495372200274L) < 1e-16);
    assert(fabsl(fabsl(m_ct[1][2]) - 4.4796575023189454e-005L) < 1e-16);

    assert(fabsl(fabsl(m_ct[2][0]) - 3.5561277024326059e-005L) < 1e-16);
    assert(fabsl(fabsl(m_ct[2][1]) - 2.8996161390351842e-005L) < 1e-16);
    assert(fabsl(fabsl(m_ct[2][2]) - 0.99999999894730907L) < 1e-16);

    return;
}

void get_celes_to_terra_matrix_test(void)
{
    double m_ct[3][3];
    get_celes_to_terra_matrix(51519.114583333336, 0.38415, m_ct);

    assert(fabs(fabs(m_ct[0][0]) - 0.44919495453022984) < 1e-15);
    assert(fabs(fabs(m_ct[0][1]) - 0.89343376515885931) < 1e-15);
    assert(fabs(fabs(m_ct[0][2]) - 9.9322034405854083e-006) < 1e-15);

    assert(fabs(fabs(m_ct[1][0]) - 0.8934337645063456) < 1e-15);
    assert(fabs(fabs(m_ct[1][1]) - 0.44919495370416457) < 1e-15);
    assert(fabs(fabs(m_ct[1][2]) - 4.4796575022978088e-005) < 1e-15);

    assert(fabs(fabs(m_ct[2][0]) - 3.5561277024326527e-005) < 1e-15);
    assert(fabs(fabs(m_ct[2][1]) - 2.8996161390348752e-005) < 1e-15);
    assert(fabs(fabs(m_ct[2][2]) - 0.99999999894730907) < 1e-15);

    return;
}


void get_terra_to_celes_matrixl_test(void)
{
    long double m_ct[3][3], m_tc[3][3];
    get_celes_to_terra_matrixl(51519.114583333336L, 0.38415L, m_ct);
    get_terra_to_celes_matrixl(m_ct, m_tc);

    assert(fabsl(fabsl(m_tc[0][0]) - 0.44919495454806802L) < 1e-16);
    assert(fabsl(fabsl(m_tc[0][1]) - 0.8934337644973771L) < 1e-16);
    assert(fabsl(fabsl(m_tc[0][2]) - 3.5561277024326059e-005L) < 1e-16);

    assert(fabsl(fabsl(m_tc[1][0]) - 0.89343376514989092L) < 1e-16);
    assert(fabsl(fabsl(m_tc[1][1]) - 0.44919495372200274L) < 1e-16);
    assert(fabsl(fabsl(m_tc[1][2]) - 2.8996161390351842e-005L) < 1e-16);

    assert(fabsl(fabsl(m_tc[2][0]) - 9.932203439718028e-006L) < 1e-16);
    assert(fabsl(fabsl(m_tc[2][1]) - 4.4796575023189454e-005L) < 1e-16);
    assert(fabsl(fabsl(m_tc[2][2]) - 0.99999999894730907L) < 1e-16);

    return;
}


void get_terra_to_celes_matrix_test(void)
{
    double m_ct[3][3], m_tc[3][3];
    get_celes_to_terra_matrix(51519.114583333336, 0.38415, m_ct);
    get_terra_to_celes_matrix(m_ct, m_tc);

    assert(fabs(fabs(m_tc[0][0]) - 0.44919495453022984) < 1e-15);
    assert(fabs(fabs(m_tc[0][1]) - 0.8934337645063456) < 1e-15);
    assert(fabs(fabs(m_tc[0][2]) - 3.5561277024326527e-005) < 1e-15);
    
    assert(fabs(fabs(m_tc[1][0]) - 0.89343376515885931) < 1e-15);
    assert(fabs(fabs(m_tc[1][1]) - 0.44919495370416457) < 1e-15);
    assert(fabs(fabs(m_tc[1][2]) - 2.8996161390348752e-005) < 1e-15);
    
    assert(fabs(fabs(m_tc[2][0]) - 9.9322034405854083e-006) < 1e-15);
    assert(fabs(fabs(m_tc[2][1]) - 4.4796575022978088e-005) < 1e-15);
    assert(fabs(fabs(m_tc[2][2]) - 0.99999999894730907) < 1e-15);

    return;
}


void spherical_to_cartesianl_test(void)
{
    long double coord[3];
    long double l = 3.787930671498L;
    long double b = -0.000003961191L;
    long double r = 148552579.338L;
    spherical_to_cartesianl(l, b, r, coord);

    assert(fabsl(fabsl(coord[0]) - 118588727.486L) < 1e-3);
    assert(fabsl(fabsl(coord[1]) - 89468332.616L) < 1e-3);
    assert(fabsl(fabsl(coord[2]) - 588.445L) < 1e-3);

    return;
}

void spherical_to_cartesian_test(void)
{
    double coord[3];
    double l = 3.787930671498;
    double b = -0.000003961191;
    double r = 148552579.338;
    spherical_to_cartesian(l, b, r, coord);

    assert(fabs(fabs(coord[0]) - 118588727.486) < 1e-3);
    assert(fabs(fabs(coord[1]) - 89468332.616) < 1e-3);
    assert(fabs(fabs(coord[2]) - 588.445) < 1e-3);

    return;
}


void test_coordinate_converters()
{
    get_mean_equ_to_fixed_matrixl_test();
    get_mean_equ_to_fixed_matrix_test();

    get_fixed_to_true_equ_matrixl_test();
    get_fixed_to_true_equ_matrix_test();

    get_true_equ_to_fixed_matrixl_test();
    get_true_equ_to_fixed_matrix_test();

    get_celes_to_terra_matrixl_test();
    get_celes_to_terra_matrix_test();

    get_terra_to_celes_matrixl_test();
    get_terra_to_celes_matrix_test();

    spherical_to_cartesianl_test();
    spherical_to_cartesian_test();
}