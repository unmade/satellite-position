//
// Created by Леша on 29.05.15.
//

#ifndef SATELLITE_PROPAGATION_SUN_H
#define SATELLITE_PROPAGATION_SUN_H

void get_sun_ecliptic_position(long double tdb, long double *l, long double *b, long double *r);


void get_sun_celestial_position(long double tdb, long double coordinates[3]);

#endif //SATELLITE_PROPAGATION_SUN_H
