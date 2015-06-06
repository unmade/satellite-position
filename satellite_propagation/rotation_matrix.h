/** \file rotation_matrix.h
 * Файл определяет матрицы поворота
 * по осям OX, OY, OZ.
 */


#ifndef SATELLITE_PROPAGATION_ROTATION_MATRIX_H
#define SATELLITE_PROPAGATION_ROTATION_MATRIX_H


/**
 * Заполняет matr матрицей поворота по оси ОX
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotxl(long double a, long double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОX
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotx(double a, double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОY
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotyl(long double a, long double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОY
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void roty(double a, double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОZ
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotzl(long double a, long double matr[3][3]);


/**
 * Заполняет matr матрицей поворота по оси ОZ
 *
 * @param[in] a Угол поворота в радианах
 * @param[out] matr Матрица поворота
 *
 */
void rotz(double a, double matr[3][3]);


/**
 * Определяет матрицу вращение Земли
 *
 * @param[in] gast Гринвическое истинное звездное время
 * @param[out] matr Матрица поворота
 *
 */
void get_earth_rotationl(long double gast, long double matr[3][3]);


/**
 * Определяет матрицу вращение Земли
 *
 * @param[in] gast Гринвическое истинное звездное время
 * @param[out] matr Матрица поворота
 *
 */
void get_earth_rotation(double gast, double matr[3][3]);


#endif //SATELLITE_PROPAGATION_ROTATION_MATRIX_H