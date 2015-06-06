/** \file propagate.h
 * Файл определяет функции для прогнозирования положения и скорости
 * небесного тела на заданный момент времени
 */

#ifndef SATELLITE_PROPAGATION_PROPAGATE_H
#define SATELLITE_PROPAGATION_PROPAGATE_H

/**
 * Выполняет прогнозирование положения небесного тела на заданную дату
 *
 * @param[in] step Шаг интегрирования
 * @param[in] start_date Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] end_date Конечный момент времени UTC, выраженный в модифицированных юлианских днях,
 *                     на который выполняется прогноз
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_vel Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 *
 */
void propagatel(long double step, long double start_date, long double end_date,
                long double start_pos[3], long double start_vel[3],
                long double final_pos[3], long double final_vel[3]);


/**
 * Выполняет прогнозирование положения небесного тела на заданную дату
 *
 * @param[in] step Шаг интегрирования
 * @param[in] start_date Начальный момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[in] end_date Конечный момент времени UTC, выраженный в модифицированных юлианских днях,
 *                     на который выполняется прогноз
 * @param[in] start_pos Начальный вектор положения в небесной СК
 * @param[in] start_vel Начальный вектор скорости в небесной СК
 * @param[out] final_pos Конечный вектор положения в небесной СК
 * @param[out] final_vel Конечный вектор скорости в небесной СК
 *
 */
void propagate(double step, double start_date, double end_date,
               double start_pos[3], double start_vel[3],
               double final_pos[3], double final_vel[3]);

#endif //SATELLITE_PROPAGATION_PROPAGATE_H
