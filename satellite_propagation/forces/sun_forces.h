//
// Created by Леша on 29.05.15.
//

#ifndef SATELLITE_PROPAGATION_SUN_FORCES_H
#define SATELLITE_PROPAGATION_SUN_FORCES_H

void get_acceleration_by_sun(long double utc_in_mjd, long double x, long double y, long double z,
                             long double acceleration[3]);

#endif //SATELLITE_PROPAGATION_SUN_FORCES_H
