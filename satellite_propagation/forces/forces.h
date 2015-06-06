/** \file forces.h
 * Файл определяет функции для вычисления ускорений,
 * обусловленных притяжением Земли, Луны, Солнца,
 * торможением атмосферы и давлением солнечного света
 */


#ifndef SATELLITE_PROPAGATION_FORCES_H
#define SATELLITE_PROPAGATION_FORCES_H

/**
 * Определяет ускорения, обусловленные притяжением действием всех сил: Земли, Луны, Солнца,
 * торможением атмосферы и давлением солнечного света
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] coord Положение тела в небесной СК
 * @param[in] vel Вектор скорости тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_forces(double utc_in_mjd, double coord[3], double vel[3], double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Земли
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_earth(double utc_in_mjd, double celes_coord[3], double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Луны
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_moon(double utc_in_mjd, double coord[3], double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Солнца
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_sun(double utc_in_mjd, double celes_coord[3], double acceleration[3]);


/**
 * Определяет ускорения, обусловленные торможением атмосферы
 *
 * @param[in] coord Вектор положения тела в небесной СК
 * @param[in] vel Вектор скорости тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_atmosphere(double coord[3], double vel[3], double acceleration[3]);


/**
 * Определяет ускорения, обусловленные солнечным давлением
 *
 * @param[in] coord Вектор положения тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_solar_pressure(double utc_in_mjd, double coord[3], double acceleration[3]);



/* Версия функций long double */

/**
 * Определяет ускорения, обусловленные притяжением действием всех сил: Земли, Луны, Солнца,
 * торможением атмосферы и давлением солнечного света
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] coord Положение тела в небесной СК
 * @param[in] vel Вектор скорости тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_forcesl(long double utc_in_mjd, long double coord[3], long double vel[3], long double acceleration[3]);


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
 * Определяет ускорения, обусловленные притяжением Луны
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_moonl(long double utc_in_mjd, long double coord[3], long double acceleration[3]);


/**
 * Определяет ускорения, обусловленные притяжением Солнца
 *
 * @param[in] utc_in_mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] celes_coord Положение тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно
 *
 */
void get_acceleration_by_sunl(long double utc_in_mjd, long double celes_coord[3],
                              long double acceleration[3]);


/**
 * Определяет ускорения, обусловленные торможением атмосферы
 *
 * @param[in] coord Вектор положения тела в небесной СК
 * @param[in] vel Вектор скорости тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_atmospherel(long double coord[3], long double vel[3], long double acceleration[3]);


/**
 * Определяет ускорения, обусловленные солнечным давлением
 *
 * @param[in] coord Вектор положения тела в небесной СК
 * @param[out] acceleration Вектор ускорений по x, y, z соответственно в небесной СК
 *
 */
void get_acceleration_by_solar_pressurel(long double utc_in_mjd, long double coord[3], long double acceleration[3]);


#endif //SATELLITE_PROPAGATION_FORCES_H
