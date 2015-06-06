//
// Created by Леша on 02.06.15.
//

#include <math.h>
#include <stdio.h>
#include "assert.h"

#include "rotation_matrix_test.h"
#include "rotation_matrix.h"

void test_rotxl(void)
{
    long double angle, matr[3][3];
    angle = 0.628318530717959L;
    rotxl(angle, matr);

    assert(matr[0][0] == 1.0L);
    assert(matr[0][1] == 0.0L);
    assert(matr[0][2] == 0.0L);

    assert(matr[1][0] == 0.0L);
    assert(matr[1][1] == cosl(angle));
    assert(matr[1][2] == sinl(angle));

    assert(matr[2][0] == 0.0L);
    assert(matr[2][1] == -sinl(angle));
    assert(matr[2][2] == cosl(angle));

//    printf("Угол a = %2.15f\n", (double)angle);
//    printf("Матрица поворота по ОХ:\n");
//    print_matrixl((long double *) matr, 3, 3);
}


void test_rotx(void)
{
    double angle, matr[3][3];
    angle = 0.628318530717959;
    rotx(angle, matr);

    assert(matr[0][0] == 1.0);
    assert(matr[0][1] == 0.0);
    assert(matr[0][2] == 0.0);

    assert(matr[1][0] == 0.0);
    assert(matr[1][1] == cos(angle));
    assert(matr[1][2] == sin(angle));

    assert(matr[2][0] == 0.0);
    assert(matr[2][1] == -sin(angle));
    assert(matr[2][2] == cos(angle));

//    printf("Угол a = %2.15f\n", angle);
//    printf("Матрица поворота по ОХ:\n");
//    print_matrix((double *) matr, 3, 3);
}


void test_rotyl(void)
{
    long double angle, matr[3][3];
    angle = 0.628318530717959L;
    rotyl(angle, matr);

    assert(matr[0][0] == cosl(angle));
    assert(matr[0][1] == 0.0L);
    assert(matr[0][2] == -sinl(angle));

    assert(matr[1][0] == 0.0L);
    assert(matr[1][1] == 1.0L);
    assert(matr[1][2] == 0.0L);

    assert(matr[2][0] == sinl(angle));
    assert(matr[2][1] == 0.0L);
    assert(matr[2][2] == cosl(angle));

//    printf("Угол a = %2.15f\n", (double)angle);
//    printf("Матрица поворота по ОY:\n");
//    print_matrixl((long double *) matr, 3, 3);
}


void test_roty(void)
{
    double angle, matr[3][3];
    angle = 0.628318530717959;
    roty(angle, matr);

    assert(matr[0][0] == cos(angle));
    assert(matr[0][1] == 0.0);
    assert(matr[0][2] == -sin(angle));

    assert(matr[1][0] == 0.0);
    assert(matr[1][1] == 1.0);
    assert(matr[1][2] == 0.0);

    assert(matr[2][0] == sin(angle));
    assert(matr[2][1] == 0.0);
    assert(matr[2][2] == cos(angle));

//    printf("Угол a = %2.15f\n", angle);
//    printf("Матрица поворота по ОY:\n");
//    print_matrix((double *) matr, 3, 3);
}


void test_rotzl(void)
{
    long double angle, matr[3][3];
    angle = 0.628318530717959L;
    rotzl(angle, matr);

    assert(matr[0][0] == cosl(angle));
    assert(matr[0][1] == sinl(angle));
    assert(matr[0][2] == 0.0L);

    assert(matr[1][0] == -sinl(angle));
    assert(matr[1][1] == cosl(angle));
    assert(matr[1][2] == 0.0L);

    assert(matr[2][0] == 0.0L);
    assert(matr[2][1] == 0.0L);
    assert(matr[2][2] == 1.0L);

//    printf("Угол a = %2.15f\n", (double)angle);
//    printf("Матрица поворота по ОZ:\n");
//    print_matrixl( (long double *) matr, 3, 3);
}


void test_rotz(void)
{
    double angle, matr[3][3];
    angle = 0.628318530717959;
    rotz(angle, matr);

    assert(matr[0][0] == cos(angle));
    assert(matr[0][1] == sin(angle));
    assert(matr[0][2] == 0.0);

    assert(matr[1][0] == -sin(angle));
    assert(matr[1][1] == cos(angle));
    assert(matr[1][2] == 0.0);

    assert(matr[2][0] == 0.0);
    assert(matr[2][1] == 0.0);
    assert(matr[2][2] == 1.0);

//    printf("Угол a = %2.15f\n", angle);
//    printf("Матрица поворота по ОZ:\n");
//    print_matrix((double *) matr, 3, 3);
}


void test_rotation_matrices(void)
{
    test_rotxl();
    test_rotx();
    test_rotyl();
    test_roty();
    test_rotzl();
    test_rotz();
}

//void test_get_earth_rotationl(void);
//void test_get_earth_rotation(void);