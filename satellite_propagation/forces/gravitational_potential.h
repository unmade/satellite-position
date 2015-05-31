/** \file gravitational_potential.h
 * Файл определяет функции предназначенные для вычисления ускорений,
 * обусловленных притяжением Земли
 */

#ifndef SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H
#define SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H

void get_acceleration_by_earth(long double utc_in_mjd, long double xc, long double yc, long double zc,
                               long double acceleration[3]);

#endif //SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H