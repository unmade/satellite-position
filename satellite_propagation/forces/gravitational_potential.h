/** \file gravitational_potential.h
 * Файл определяет функции предназначенные для вычисления ускорений,
 * обусловленных притяжением Земли
 */

#ifndef SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H
#define SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H

/**
 * Определяет ускорения, обусловленные притяжением Земли
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_earthl(long double utc_in_mjd, long double celes_coord[3],
                                long double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Земли
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_earth(double utc_in_mjd, double celes_coord[3], double acceleration[3]);

#endif //SATELLITE_PROPAGATION_GRAVITATIONAL_POTENTIAL_H