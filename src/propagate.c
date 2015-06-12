//
// Created by Леша on 31.05.15.
//

#include <math.h>
#include "integration.h"
#include "propagate.h"

#include "constants.h"

void propagatel(int method, long double step, long double start_date, long double end_date,
                long double start_pos[3], long double start_vel[3],
                long double final_pos[3], long double final_vel[3])
{
    int i, j;
    long double a[7][3];
    long double pos[3], vel[3];
    long double date_diff = (end_date - start_date) * SEC_IN_DAY;
    long double count = truncl(date_diff / step);
    long double fin_step = date_diff - count*step;
    long double add = step / SEC_IN_DAY;
    
    pos[0] = start_pos[0];
    pos[1] = start_pos[1];
    pos[2] = start_pos[2];
    vel[0] = start_vel[0];
    vel[1] = start_vel[1];
    vel[2] = start_vel[2];

    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            a[i][j] = 0;

    i = 0;
    
    while (i < count)
    {
        if (method == 0)
            everhartl(start_date, step, pos, vel, final_pos, final_vel, a);
        else if (method == 1)
            rungekuttal(start_date, step, pos, vel, final_pos, final_vel);
        else
            return;
        pos[0] = final_pos[0];
        pos[1] = final_pos[1];
        pos[2] = final_pos[2];
        vel[0] = final_vel[0];
        vel[1] = final_vel[1];
        vel[2] = final_vel[2];
        start_date += add;
        ++i;
    }
    start_date -= add;
    if (fin_step)
    {
        start_date += fin_step / SEC_IN_DAY;
        if (method == 0)
            everhartl(start_date, fin_step, pos, vel, final_pos, final_vel, a);
        else if (method == 1)
            rungekuttal(start_date, fin_step, pos, vel, final_pos, final_vel);
        else
            return;
    }

    return;
}


void propagate(int method, double step, double start_date, double end_date,
               double start_pos[3], double start_vel[3],
               double final_pos[3], double final_vel[3])
{
    int i, j;
    double a[7][3];
    double pos[3], vel[3];
    double date_diff = (end_date - start_date) * SEC_IN_DAY;
    double count = trunc(date_diff / step);
    double fin_step = date_diff - count*step;
    double add = step / SEC_IN_DAY;

    pos[0] = start_pos[0];
    pos[1] = start_pos[1];
    pos[2] = start_pos[2];
    vel[0] = start_vel[0];
    vel[1] = start_vel[1];
    vel[2] = start_vel[2];

    for (i = 0; i < 7; i++)
        for (j = 0; j < 3; j++)
            a[i][j] = 0;

    i = 0;
    while (i < count)
    {
        if (method == 0)
            everhart(start_date, step, pos, vel, final_pos, final_vel, a);
        else if (method == 1)
            rungekutta(start_date, step, pos, vel, final_pos, final_vel);
        else
            return;
        pos[0] = final_pos[0];
        pos[1] = final_pos[1];
        pos[2] = final_pos[2];
        vel[0] = final_vel[0];
        vel[1] = final_vel[1];
        vel[2] = final_vel[2];
        start_date += add;
        ++i;
    }
    start_date -= add;
    if (fin_step)
    {
        start_date += fin_step / SEC_IN_DAY;
        if (method == 0)
            everhart(start_date, fin_step, pos, vel, final_pos, final_vel, a);
        else if (method == 1)
            rungekutta(start_date, fin_step, pos, vel, final_pos, final_vel);
        else
            return;
    }

    return;
}