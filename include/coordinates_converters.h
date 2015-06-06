/** \file coordinates_converters.h
 * Файл содержит функции для перехода между СК
 */


#ifndef SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H
#define SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H

/**
 * Функция реализует переход от средней подвижной экваториальной СК к
 * небесной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[out] m_p[3][3] Матрица преобразования от средней подвижной экваториальной СК
 * к небесной
 *
 */
void get_mean_equ_to_fixed_matrixl(long double precession_matrix[3][3], long double m_p[3][3]);


/**
 * Функция реализует переход от средней подвижной экваториальной СК к
 * небесной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[out] m_p[3][3] Матрица преобразования от средней подвижной экваториальной СК
 * к небесной
 *
 */
void get_mean_equ_to_fixed_matrix(double precession_matrix[3][3], double m_p[3][3]);


/**
 * Функция реализует переход от небесной к истинной экваториальной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[in] nutation_matrix[3][3] Матрица нутации
 * @param[out] m_np[3][3] Матрица преобразования от небесной к истинной экваториальной СК
 *
 */
void get_fixed_to_true_equ_matrixl(long double precession_matrix[3][3], long double nutation_matrix[3][3],
                                   long double m_np[3][3]);


/**
 * Функция реализует переход от небесной к истинной экваториальной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[in] nutation_matrix[3][3] Матрица нутации
 * @param[out] m_np[3][3] Матрица преобразования от небесной к истинной экваториальной СК
 *
 */
void get_fixed_to_true_equ_matrix(double precession_matrix[3][3], double nutation_matrix[3][3],
                                  double m_np[3][3]);


/**
 * Функция реализует переход от истинной экваториальной СК к небесной
 *
 * @param[in] m_np Матрица преобразования от небесной к истинной экваториальной СК
 * @param[out] m_pnМатрица преобразования от истинной экваториальной СК к небесной
 *
 */
void get_true_equ_to_fixed_matrixl(long double m_np[3][3], long double m_pn[3][3]);


/**
 * Функция реализует переход от истинной экваториальной СК к небесной
 *
 * @param[in] m_np Матрица преобразования от небесной к истинной экваториальной СК
 * @param[out] m_pnМатрица преобразования от истинной экваториальной СК к небесной
 *
 */
void get_true_equ_to_fixed_matrix(double m_np[3][3], double m_pn[3][3]);


/**
 * Функция реализует переход от небесной СК в земную
 *
 * @param[in] utc_in_mjd Время UTC в модифицированных юлианских днях
 * @param[in] delta_ut Поправка ∆T = UT1 - UTC
 * @param[out] m_ct Матрица перехода от небесной к земной СК
 *
 */
void get_celes_to_terra_matrixl(long double utc_in_mjd, long double delta_ut, long double m_ct[3][3]);


/**
 * Функция реализует переход от небесной СК в земную
 *
 * @param[in] utc_in_mjd Время UTC в модифицированных юлианских днях
 * @param[in] delta_ut Поправка ∆T = UT1 - UTC
 * @param[out] m_ct Матрица перехода от небесной к земной СК
 *
 */
void get_celes_to_terra_matrix(double utc_in_mjd, double delta_ut, double m_ct[3][3]);


/**
 * Функция реализует переход из земной СК в небесную
 *
 * @param[in] m_ct Матрица перехода от небесной к земной СК
 * @param[out] m_tc Матрица преобразования от истинной экваториальной СК к небесной
 *
 */
void get_terra_to_celes_matrixl(long double m_ct[3][3], long double m_tc[3][3]);


/**
 * Функция реализует переход из земной СК в небесную
 *
 * @param[in] m_ct Матрица перехода от небесной к земной СК
 * @param[out] m_tc Матрица преобразования от истинной экваториальной СК к небесной
 *
 */
void get_terra_to_celes_matrix(double m_ct[3][3], double m_tc[3][3]);


/**
 * Функция реализует переход из сферической (эклиптической) в декартову (прямоугольную) СК
 *
 * @param[in] l Долгота [рад]
 * @param[in] b Широта  [рад]
 * @param[in] r Расстояние  [километр]
 * @param[out] coordinates Содержит прямоугольные коор-ты x, y, z соответственно [километр]
 *
 */
void spherical_to_cartesianl(long double l, long double b, long double r,
                             long double coordinates[3]);


/**
 * Функция реализует переход из сферической (эклиптической) в декартову (прямоугольную) СК
 *
 * @param[in] l Долгота [рад]
 * @param[in] b Широта  [рад]
 * @param[in] r Расстояние  [километр]
 * @param[out] coordinates Содержит прямоугольные коор-ты x, y, z соответственно [километр]
 *
 */
void spherical_to_cartesian(double l, double b, double r, double coordinates[3]);


/**
 * Функция реализует переход из геодезической в декартову (прямоугольную) СК
 *
 * @param[in] phi Широта [рад]
 * @param[in] lambda Долгота  [рад]
 * @param[in] h Высота  [километр]
 * @param[out] coordinates Содержит прямоугольные коор-ты x, y, z соответственно [километр]
 *
 */
void geodesic_to_cartesianl(long double phi, long double lambda, long double h,
                            long double coordinates[3]);


/**
 * Функция реализует переход из геодезической в декартову (прямоугольную) СК
 *
 * @param[in] phi Широта [рад]
 * @param[in] lambda Долгота  [рад]
 * @param[in] h Высота  [километр]
 * @param[out] coordinates Содержит прямоугольные коор-ты x, y, z соответственно [километр]
 *
 */
void geodesic_to_cartesian(double phi, double lambda, double h,
                           double coordinates[3]);

#endif //SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H