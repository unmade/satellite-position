//
// Created by Леша on 22.05.15.
//
#include "math.h"
#include <stdlib.h>

#include "rotation_matrix.h"


void rotate_by_x(long double a, long double matr[3][3])
{
    matr[0][0] = 1;
    matr[0][1] = 0;
    matr[0][2] = 0;

    matr[1][0] = 0;
    matr[1][1] = cosl(a);
    matr[1][2] = sinl(a);

    matr[2][0] = 0;
    matr[2][1] = -sinl(a);
    matr[2][2] = cosl(a);

    return;
}


void rotate_by_y(long double a, long double matr[3][3])
{
    matr[0][0] = cosl(a);
    matr[0][1] = 0;
    matr[0][2] = -sinl(a);

    matr[1][0] = 0;
    matr[1][1] = 1;
    matr[1][2] = 0;

    matr[2][0] = sinl(a);
    matr[2][1] = 0;
    matr[2][2] = cosl(a);

    return;
}


void rotate_by_z(long double a, long double matr[3][3])
{
    matr[0][0] = cosl(a);
    matr[0][1] = sinl(a);
    matr[0][2] = 0;

    matr[1][0] = -sinl(a);
    matr[1][1] = cosl(a);
    matr[1][2] = 0;

    matr[2][0] = 0;
    matr[2][1] = 0;
    matr[2][2] = 1;

    return;
}


void get_earth_rotation_matrix(long double gast, long double matr[3][3])
{
    rotate_by_z(gast, matr);
    return;
}


//void destroyArray(long double** arr)
//{
//    free(*arr);
//    free(arr);
//}


//long double **rotate_by_y(long double a)
//{
//    int i;
//    long double **matr = (long double **)malloc(N * sizeof(long double *));
//    for (i=0; i<N; i++)
//        matr[i] = (long double *)malloc(N * sizeof(long double));
//
//    matr[0][0] = cos(a);
//    matr[0][1] = 0;
//    matr[0][2] = -sin(a);
//
//    matr[1][0] = 0;
//    matr[1][1] = 1;
//    matr[1][2] = 0;
//
//    matr[2][0] = sin(a);
//    matr[2][1] = 0;
//    matr[2][2] = cos(a);
//
//    return matr;

//    long double matr[N][N] = { {cos(a), 0, -sin(a)},
//                          {0, 1, 0},
//                          {sin(a), 0, cos(a)}};
//}