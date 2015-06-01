//
// Created by Леша on 30.05.15.
//

#ifndef SATELLITE_PROPAGATION_EVERHART_H
#define SATELLITE_PROPAGATION_EVERHART_H

void do_everhart(long double utc_in_mjd, long double start_pos[3], long double final_pos[3], long double start_vel[3],
                 long double fin_vel[3], long double a[7][3], long double step);

#endif //SATELLITE_PROPAGATION_EVERHART_H
