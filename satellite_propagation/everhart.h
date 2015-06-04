//
// Created by Леша on 30.05.15.
//

#ifndef SATELLITE_PROPAGATION_EVERHART_H
#define SATELLITE_PROPAGATION_EVERHART_H

/**
 * Функция реализует метод Эверхарта (неявный одношаговый метод)
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
void do_everhartl(long double utc_in_mjd, long double step, long double start_pos[3], long double start_vel[3],
                  long double final_pos[3], long double fin_vel[3], long double alpha[7][3]);


/**
 * Функция реализует метод Эверхарта (неявный одношаговый метод)
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
void do_everhart(double utc_in_mjd, double step, double start_pos[3], double start_vel[3], double final_pos[3 ],
                 double fin_vel[3], double alpha[7][3]);

#endif //SATELLITE_PROPAGATION_EVERHART_H
