#include <stdio.h>
#include "matrix_operations.h"


void mult_matricesl(long double a[3][3], long double b[3][3], long double c[3][3])
{
    int i, j, k;
    for (i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            c[i][j] = 0;
            for (k=0; k<3; k++)
            {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}


void mult_matrices(double a[3][3], double b[3][3], double c[3][3])
{
    int i, j, k;
    for (i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            c[i][j] = 0;
            for (k=0; k<3; k++)
            {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}


void mult_matrix_by_vectorl(long double matr[3][3], long double vec[3], long double res[3])
{
    int i, j;
    for (i=0; i<3; i++)
    {
        res[i] = 0;
        for (j=0; j<3; j++)
        {
            res[i] += matr[i][j]*vec[j];
        }
    }
}



void mult_matrix_by_vector(double matr[3][3], double vec[3], double res[3])
{
    int i, j;
    for (i=0; i<3; i++)
    {
        res[i] = 0;
        for (j=0; j<3; j++)
        {
            res[i] += matr[i][j]*vec[j];
        }
    }
}



void transposel(void *dest, void *src, int src_h, int src_w)
{
    int i, j;
    long double (*d)[src_h] = dest, (*s)[src_w] = src;
    for (i = 0; i < src_h; i++)
        for (j = 0; j < src_w; j++)
            d[j][i] = s[i][j];
}


void transpose(void *dest, void *src, int src_h, int src_w)
{
    int i, j;
    double (*d)[src_h] = dest, (*s)[src_w] = src;
    for (i = 0; i < src_h; i++)
        for (j = 0; j < src_w; j++)
            d[j][i] = s[i][j];
}


void print_matrixl(long double *matrix, int n, int m)
{
    int i, j;
    for (i=0; i<n; i++)
    {
        for (j=0;j<m;j++)
        {
            printf("%2.15f   ", (double)(*((matrix+i*n) + j)));
        }
        printf("\n");
    }
}


void print_matrix(double *matrix, int n, int m)
{
    int i, j;
    for (i=0; i<n; i++)
    {
        for (j=0;j<m;j++)
        {
            printf("%2.15f   ", *((matrix+i*n) + j));
        }
        printf("\n");
    }
}