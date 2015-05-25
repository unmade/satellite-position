/** \file rotation_matrix.h
 * Файл определяет матрицы поворота
 * по осям OX, OY, OZ.
 */


#ifndef SATELLITE_PROPAGATION_ROTATION_MATRIX_H
#define SATELLITE_PROPAGATION_ROTATION_MATRIX_H
#endif //SATELLITE_PROPAGATION_ROTATION_MATRIX_H


/**
 * Заполняет matr матрицей поворота по оси ОX
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotate_by_x(double a, double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОY
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotate_by_y(double a, double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОZ
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotate_by_z(double a, double matr[3][3]);


/**
 * Определяет матрицу вращение Земли
 *
 * @param[in] gast Гринвическое истинное звездное время
 * @param[out] matr Матрица поворота
 *
 */
void get_earth_rotation_matrix(double gast, double matr[3][3]);


//double** rotate_by_y(double a);

//void destroyArray(double** arr);
