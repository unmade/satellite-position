//
// Created by Леша on 22.05.15.
//
#include "math.h"

#include "rotation_matrix.h"


void rotxl(long double a, long double matr[3][3])
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


void rotx(double a, double matr[3][3])
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


void rotyl(long double a, long double matr[3][3])
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


void roty(double a, double matr[3][3])
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


void rotzl(long double a, long double matr[3][3])
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


void rotz(double a, double matr[3][3])
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


void get_earth_rotationl(long double gast, long double matr[3][3])
{
    rotzl(gast, matr);
    return;
}


void get_earth_rotation(double gast, double matr[3][3])
{
    rotz(gast, matr);
    return;
}
