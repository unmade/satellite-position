//
// Created by Леша on 22.05.15.
//
#include "math.h"
#include <stdlib.h>

#include "rotation_matrix.h"


int rotate_by_x(double a, double matr[N][N])
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

    return 0;
}


int rotate_by_y(double a, double matr[N][N]) {
    matr[0][0] = cos(a);
    matr[0][1] = 0;
    matr[0][2] = -sin(a);

    matr[1][0] = 0;
    matr[1][1] = 1;
    matr[1][2] = 0;

    matr[2][0] = sin(a);
    matr[2][1] = 0;
    matr[2][2] = cos(a);

    return 0;
}

int rotate_by_z(double a, double matr[N][N])
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

    return 0;
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