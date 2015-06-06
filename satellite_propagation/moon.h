/** \file moon.h
 * Файл определяет функции для вычисления положений Луны
 */

#ifndef SATELLITE_PROPAGATION_MOON_H
#define SATELLITE_PROPAGATION_MOON_H



/**
 * Функция вычисляет геоцентрические эклиптические координаты Луны
 *
 * @param[in] tdb Барицентрическое динамическое время в модифицированных юлианских днях
 * @param[out] l долгота
 * @param[out] b широта
 * @param[out] r расстояние
 *
 */
void get_moon_ecliptic_positionl(long double tdb, long double *l, long double *b, long double *r);


/**
 * Функция вычисляет геоцентрические эклиптические координаты Луны
 *
 * @param[in] tdb Барицентрическое динамическое время в модифицированных юлианских днях
 * @param[out] l Долгота
 * @param[out] b Широта
 * @param[out] r Расстояние
 *
 */
void get_moon_ecliptic_position(double tdb, double *l, double *b, double *r);


/**
 * Функция вычисляет положение Луны в небесной системе координат
 *
 * @param[in] tdb Барицентрическое динамическое время в модифицированных юлианских днях
 * @param[out] coordinates Прямоугольные координаты x, y, z соответственно
 *
 */
void get_moon_celestial_positionl(long double tdb, long double coordinates[3]);


/**
 * Функция вычисляет положение Луны в небесной системе координат
 *
 * @param[in] tdb Барицентрическое динамическое время в модифицированных юлианских днях
 * @param[out] coordinates Прямоугольные координаты x, y, z соответственно
 *
 */
void get_moon_celestial_position(double tdb, double coordinates[3]);

#endif //SATELLITE_PROPAGATION_MOON_H