//
// Created by Леша on 31.05.15.
//

#include <math.h>
#include "propagate.h"
#include "everhart.h"
#include "constants.h"

void propagatel(int step, long double start_date, long double end_date,
                long double start_pos[3], long double start_vel[3],
                long double final_pos[3], long double final_vel[3])
{
    long double date_diff = (end_date - start_date) * SEC_IN_DAY;
    long double count = truncl(date_diff / step);
    long double fin_step = date_diff - count*step;

    int i, j;

    long double a[7][3];
    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            a[i][j] = 0;

    long double add = step / SEC_IN_DAY;

    i = 0;
    while (i < count)
    {
        do_everhartl(start_date, step, start_pos, start_vel, final_pos, final_vel, a);

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
        do_everhartl(start_date, fin_step, start_pos, start_vel, final_pos, final_vel, a);
    }

    return;
}


void propagate(int step, double start_date, double end_date,
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
        do_everhart(start_date, step, start_pos, start_vel, final_pos, final_vel, a);

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
        do_everhart(start_date, fin_step, start_pos, start_vel, final_pos, final_vel, a);
    }

    return;
}