//
// Created by Леша on 06.06.15.
//

#ifndef SATELLITE_PROPAGATION_INTEGRATION_H
#define SATELLITE_PROPAGATION_INTEGRATION_H

/**
 * Функция реализует метод Эверхарта 15-го порядка (неявный одношаговый метод)
 *
 * @param[in] utc_in_mjd Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] step Шаг интегрирования
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_pos Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 * @param[in,out] alpha Коэффициенты альфа
 *
 */
void everhart(double utc_in_mjd, double step, double start_pos[3], double start_vel[3], double final_pos[3],
              double fin_vel[3], double alpha[7][3]);


/**
 * Функция реализует метод Рунге-Кутты 4-го порядка
 *
 * @param[in] utc_in_mjd Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] step Шаг интегрирования
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_vel Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 *
 */
void rungekutta(double mjd, double step,
                double r[3], double v[3],
                double fin_r[3], double fin_v[3]);


/**
 * Функция реализует метод Эверхарта 15-го порядка (неявный одношаговый метод)
 *
 * @param[in] utc_in_mjd Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] step Шаг интегрирования
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_vel Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 * @param[in,out] alpha Коэффициенты альфа
 *
 */
void everhartl(long double utc_in_mjd, long double step, long double start_pos[3], long double start_vel[3],
               long double final_pos[3], long double fin_vel[3], long double alpha[7][3]);


/**
 * Функция реализует метод Рунге-Кутты 4-го порядка
 *
 * @param[in] utc_in_mjd Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] step Шаг интегрирования
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_vel Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 *
 */
void rungekuttal(long double mjd, long double step,
                 long double r[3], long double v[3],
                 long double fin_r[3], long double fin_v[3]);


#endif //SATELLITE_PROPAGATION_INTEGRATION_H
