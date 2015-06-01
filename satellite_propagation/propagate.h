//
// Created by Леша on 31.05.15.
//

#ifndef SATELLITE_PROPAGATION_PROPAGATE_H
#define SATELLITE_PROPAGATION_PROPAGATE_H

void propagate(int step, long double start_date, long double end_date,
               long double start_pos[3], long double start_vel[3],
               long double final_pos[3], long double end_vel[3]);

#endif //SATELLITE_PROPAGATION_PROPAGATE_H
