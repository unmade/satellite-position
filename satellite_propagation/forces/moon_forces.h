//
// Created by Леша on 29.05.15.
//

#ifndef SATELLITE_PROPAGATION_MOON_FORCES_H
#define SATELLITE_PROPAGATION_MOON_FORCES_H

/**
 * Определяет ускорения, обусловленные притяжением Луны
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_moonl(long double utc_in_mjd, long double coord[3], long double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Луны
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_moon(double utc_in_mjd, double coord[3], double acceleration[3]);

#endif //SATELLITE_PROPAGATION_MOON_FORCES_H
