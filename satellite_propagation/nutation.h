/** \file nutation.h
 * Файл определяет функции предназначенные для вычисления
 * матрицы нутации и её параметров
 */
#ifndef SATELLITE_PROPAGATION_NUTATION_MATRIX_H
#define SATELLITE_PROPAGATION_NUTATION_MATRIX_H

#endif //SATELLITE_PROPAGATION_NUTATION_MATRIX_H

/**
 * Вычисляет числовое значение ε(t) угла наклона эклиптики к среднему экватору
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени
 *
 *  @return Числовое значение ε(t)
 */
long double get_eps_mean(long double tdb);


/**
 * Вычисляет числовые значения фундаментальных аргументов
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени
 * @param[out] fund_args[5] Массив с фундаментальными аргументами.
 * Фундаментальные аргументы идут в следующем порядке: λ(t), l(t), l'(t), F(t), D(t)
 */
void get_fundumental_args(long double tdb, long double fund_args[5]);


/**
 * Вычисляет скорректированные числовые значения фундаментальных аргументов
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени
 * @param[out] corr_fund_args[5] Массив с фундаментальными аргументами.
 * Фундаментальные аргументы идут в следующем порядке: λ(t), l(t), l'(t), F(t), D(t)
 */
void get_corr_fundumental_args(long double tdb, long double corr_fund_args[5]);


/**
 * Вычисляет параметры нутации
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени
 * @param[out] phi Нутация в долготе
 * @param[out] eps Нутация в наклоне
 *
 */
void get_nutation_parameters(long double tdb, long double *psi, long double *eps);


/**
 * Вычисляет матрицу нутации
 *
 * @param[in] tdb Момент времени в шкале барицентрического динамического времени
 * @param[out] nutation_matrix[3][3] Матрица нутации
 *
 */
void get_nutation_matrix(long double tdb, long double nutation_matrix[3][3]);