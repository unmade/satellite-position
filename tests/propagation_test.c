//
// Created by user on 04.06.2015.
//

#include <math.h>
#include <assert.h>
#include <integration.h>
#include "propagation_test.h"
#include "date_converters.h"
#include "propagate.h"

void propagatel_test(void)
{
    long double start_date = utc_to_mjdl(2015, 5, 14, 6, 0, 0);
    long double end_date = utc_to_mjdl(2015, 5, 15, 6, 0, 0.0);
    long double pos[3], fin_pos[3], vel[3], fin_vel[3];

    pos[0] = 33508.4071859207L;
    pos[1] = 25585.7790086467L;
    pos[2] = -493.907358861003L;
    vel[0] = -1.86508722606082L;
    vel[1] = 2.44434537660644L;
    vel[2] = 0.0378635308063818L;

    propagatel(0, 62.0L, start_date, end_date, pos, vel, fin_pos, fin_vel);

    assert(fabsl(fabsl(fin_pos[0]) - 33066.259498924905L) < 1e-4);
    assert(fabsl(fabsl(fin_pos[1]) - 26154.874481927356L) < 1e-4);
    assert(fabsl(fabsl(fin_pos[2]) - 484.72224379488227L) < 1e-4);

    return;
}

void propagate_test(void)
{
    double start_date = utc_to_mjd(2015, 5, 14, 6, 0, 0);
    double end_date = utc_to_mjd(2015, 5, 15, 6, 0, 0.0);
    double pos[3], fin_pos[3], vel[3], fin_vel[3];

    pos[0] = 33508.4071859207;
    pos[1] = 25585.7790086467;
    pos[2] = -493.907358861003;
    vel[0] = -1.86508722606082;
    vel[1] = 2.44434537660644;
    vel[2] = 0.0378635308063818;

    propagate(0, 62.0, start_date, end_date, pos, vel, fin_pos, fin_vel);

    assert(fabs(fabs(fin_pos[0]) - 33066.259498923195) < 1e-4);
    assert(fabs(fabs(fin_pos[1]) - 26154.874481928142) < 1e-4);
    assert(fabs(fabs(fin_pos[2]) - 484.72224379478911) < 1e-4);

    return;
}


void test_propagation(void)
{
    propagatel_test();
    propagate_test();
    return;
}