//
// Created by Леша on 31.05.15.
//

#include <math.h>
#include "propagate.h"
#include "integration/integration.h"
#include "constants.h"

void propagatel(long double step, long double start_date, long double end_date,
                long double start_pos[3], long double start_vel[3],
                long double final_pos[3], long double final_vel[3])
{
    int i, j;
    long double a[7][3];
    long double pos[3], vel[3];
    long double date_diff = (end_date - start_date) * SEC_IN_DAY;
    long double count = truncl(date_diff / step);
    long double fin_step = date_diff - count*step;

    pos[0] = start_pos[0];
    pos[1] = start_pos[1];
    pos[2] = start_pos[2];
    vel[0] = start_vel[0];
    vel[1] = start_vel[1];
    vel[2] = start_vel[2];

    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            a[i][j] = 0;

    long double add = step / SEC_IN_DAY;

    i = 0;
    while (i < count)
    {
        everhartl(start_date, step, pos, vel, final_pos, final_vel, a);
//        rungekuttal(start_date, step, pos, vel, final_pos, final_vel);

        pos[0] = final_pos[0];
        pos[1] = final_pos[1];
        pos[2] = final_pos[2];
        vel[0] = final_vel[0];
        vel[1] = final_vel[1];
        vel[2] = final_vel[2];
        start_date += add;
        ++i;
    }

    if (fin_step)
    {
        start_date += fin_step / SEC_IN_DAY;
//        rungekuttal(start_date, fin_step, pos, vel, final_pos, final_vel);
        everhartl(start_date, fin_step, pos, vel, final_pos, final_vel, a);
    }

    return;
}


void propagate(double step, double start_date, double end_date,
               double start_pos[3], double start_vel[3],
               double final_pos[3], double final_vel[3])
{
    double date_diff = (end_date - start_date) * SEC_IN_DAY;
    double count = trunc(date_diff / step);
    double fin_step = date_diff - count*step;

    int i, j;

    double a[7][3];
    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            a[i][j] = 0;

    double add = step / SEC_IN_DAY;

    i = 0;
    while (i < count)
    {
        everhart(start_date, step, start_pos, start_vel, final_pos, final_vel, a);


        start_pos[0] = final_pos[0];
        start_pos[1] = final_pos[1];
        start_pos[2] = final_pos[2];
        start_vel[0] = final_vel[0];
        start_vel[1] = final_vel[1];
        start_vel[2] = final_vel[2];
        start_date += add;
        ++i;
    }

    if (fin_step)
    {
        start_date += fin_step / SEC_IN_DAY;
        everhart(start_date, fin_step, start_pos, start_vel, final_pos, final_vel, a);
    }

    return;
}