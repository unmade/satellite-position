/** \file precession_matrix.h
 * Файл содержит функции для вычисления матрицы прецессии
 * и её параметров.
 */
#ifndef SATELLITE_PROPAGATION_PRECESSION_MATRIX_H
#define SATELLITE_PROPAGATION_PRECESSION_MATRIX_H

#endif //SATELLITE_PROPAGATION_PRECESSION_MATRIX_H

/**
 * Вычисляет параметры прецессии
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени TDB,
 * выраженный в модифицированных юлианских днях
 * @param[out] zeta Угловая переменная ζ(t)
 * @param[out] theta Угловая переменная θ(t)
 * @param[out] z Угловая переменная z(t)
 *
 */
void calc_precession_parameters(double tdb, double *zeta, double *theta, double *z);


/**
 * Вычисляет матрицу прецессии
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени TDB,
 * выраженный в модифицированных юлианских днях
 * @param[out] zeta Угловая переменная ζ(t)
 * @param[out] theta Угловая переменная θ(t)
 * @param[out] z Угловая переменная z(t)
 *
 */
void calc_precession_matrix(double zeta, double theta, double z,
                            double precession_matr[3][3]);
