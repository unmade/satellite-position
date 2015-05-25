//
// Created by Леша on 22.05.15.
//
#include "math.h"
#include <stdlib.h>

#include "rotation_matrix.h"


void rotate_by_x(double a, double matr[3][3])
{
    matr[0][0] = 1;
    matr[0][1] = 0;
    matr[0][2] = 0;

    matr[1][0] = 0;
    matr[1][1] = cos(a);
    matr[1][2] = sin(a);

    matr[2][0] = 0;
    matr[2][1] = -sin(a);
    matr[2][2] = cos(a);

    return;
}


void rotate_by_y(double a, double matr[3][3])
{
    matr[0][0] = cos(a);
    matr[0][1] = 0;
    matr[0][2] = -sin(a);

    matr[1][0] = 0;
    matr[1][1] = 1;
    matr[1][2] = 0;

    matr[2][0] = sin(a);
    matr[2][1] = 0;
    matr[2][2] = cos(a);

    return;
}


void rotate_by_z(double a, double matr[3][3])
{
    matr[0][0] = cos(a);
    matr[0][1] = sin(a);
    matr[0][2] = 0;

    matr[1][0] = -sin(a);
    matr[1][1] = cos(a);
    matr[1][2] = 0;

    matr[2][0] = 0;
    matr[2][1] = 0;
    matr[2][2] = 1;

    return;
}


void get_earth_rotation_matrix(double gast, double matr[3][3])
{
    rotate_by_z(gast, matr);
    return;
}


//void destroyArray(double** arr)
//{
//    free(*arr);
//    free(arr);
//}


//double **rotate_by_y(double a)
//{
//    int i;
//    double **matr = (double **)malloc(N * sizeof(double *));
//    for (i=0; i<N; i++)
//        matr[i] = (double *)malloc(N * sizeof(double));
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

//    double matr[N][N] = { {cos(a), 0, -sin(a)},
//                          {0, 1, 0},
//                          {sin(a), 0, cos(a)}};
//}