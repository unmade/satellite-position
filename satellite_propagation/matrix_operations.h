/** \file matrix_operations.h
 * Файл определяет функции реализующие операции над матрицами
 */

#ifndef SATELLITE_PROPAGATION_MATRIX_OPERATIONS_H
#define SATELLITE_PROPAGATION_MATRIX_OPERATIONS_H

#endif //SATELLITE_PROPAGATION_MATRIX_OPERATIONS_H


/**
 * Выполняет операцию произведения 2-х матриц размером 3x3
 *
 * @param[in] a Первая матрица
 * @param[in] b Вторая матрица
 * @param[out] c Матрица, содержащая результат произведения a*b
 *
 */
void mult_matrix(long double a[3][3], long double b[3][3], long double c[3][3]);


/**
 * Выполняет операцию произведения матрицы на вектор
 *
 * @param[in] a Матрица
 * @param[in] vec Вектор
 * @param[out] res Результат произведения (вектор)
 *
 */
void mult_matrix_by_vector(long double matr[3][3], long double vec[3], long double res[3]);


/**
 * Выполняет операцию транспонирования матрицы
 *
 * @param[out] *dest Транспонированная матрица
 * @param[int] *src Исходная матрица
 * @param[in] src_h Количество столбцов в матрице
 * @param[in] src_w Количество строк в матрице
 *
 */
void transpose(void *dest, void *src, int src_h, int src_w);

/**
 * Выводит матрицу на печать
 *
 * @param[in] n Количество строк в матрицу
 * @param[in] m Количество столбцов в матрице
 * @param[in] *matrix Указатель на матрицу, которую нужно вывести на печать
 *
 */
void print_matrix(int n, int m, long double *matrix);