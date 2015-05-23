/** \file rotation_matrix.h
 * Файл определяет матрицы поворота
 * по осям OX, OY, OZ.
 */


#ifndef SATELLITE_PROPAGATION_ROTATION_MATRIX_H
#define SATELLITE_PROPAGATION_ROTATION_MATRIX_H
#define N 3
#endif //SATELLITE_PROPAGATION_ROTATION_MATRIX_H



/**
 * Заполняет matr матрицей поворота по оси ОX
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 * @return 0
 */
int rotate_by_x(double a, double matr[N][N]);


/**
 * Заполняет matr матрицей поворота по оси ОY
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 * @return 0
 */
int rotate_by_y(double a, double matr[N][N]);


/**
 * Заполняет matr матрицей поворота по оси ОZ
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 * @return 0
 */
int rotate_by_z(double a, double matr[N][N]);


//double** rotate_by_y(double a);

//void destroyArray(double** arr);
