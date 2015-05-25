/** \file coordinates_converters.h
 * Файл содержит функции для перехода между СК
 */


#ifndef SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H
#define SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H

#endif //SATELLITE_PROPAGATION_COORDINATES_CONVERTERS_H

/**
 * Функция реализует переход от средней подвижной экваториальной СК к
 * небесной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[out] m_p[3][3] Матрица преобразования от средней подвижной экваториальной СК
 * к небесной
 *
 */
void mean_equ_to_fixed(double precession_matrix[3][3], double m_p[3][3]);


/**
 * Функция реализует переход от небесной к истинной экваториальной СК
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[in] nutation_matrix[3][3] Матрица нутации
 * @param[out] m_np[3][3] Матрица преобразования от небесной к истинной экваториальной СК
 *
 */
void fixed_to_true_equ(double precession_matrix[3][3], double nutation_matrix[3][3],
                       double m_np[3][3]);


/**
 * Функция реализует переход от истинной экваториальной СК к небесной
 *
 * @param[in] m_np[3][3] Матрица преобразования от небесной к истинной экваториальной СК
 * @param[out] m_pn[3][3] Матрица преобразования от истинной экваториальной СК к небесной
 *
 */
void true_equ_to_fixed(double m_np[3][3], double m_pn[3][3]);

/**
 * Функция реализует переход от небесной СК в земную
 *
 * @param[in] precession_matrix[3][3] Матрица прецессии
 * @param[in] nutation_matrix[3][3] Матрица нутации
 * @param[in] earth_rotation_matrix Матрица вращения Земли
 * @param[out] m_ct[3][3] Матрица перехода от небесной к земной СК
 *
 */
void fixed_to_terra(double precession_matrix[3][3], double nutation_matrix[3][3],
                    double earth_rotation_matrix[3][3], double m_ct[3][3]);


/**
 * Функция реализует переход из земной СК в небесную
 *
 * @param[in] m_ct[3][3] Матрица перехода от небесной к земной СК
 * @param[out] m_tc[3][3] Матрица преобразования от истинной экваториальной СК к небесной
 *
 */
void terra_to_fixed(double m_ct[3][3], double m_tc[3][3]);