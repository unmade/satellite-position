/** \file gravitational_potential.h
 * Файл определяет функции предназначенные для вычисления ускорений,
 * обусловленных притяжением Земли
 */

#ifndef SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H
#define SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H

#endif //SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H

void calc(double xc, double yc, double zc,
          double utc_in_mjd,
          double *Fx, double *Fy, double *Fz);