#include <stdio.h>
#include "matrix_operations.h"



void mult_matrix(double a[3][3], double b[3][3], double c[3][3])
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


void transpose(void *dest, void *src, int src_h, int src_w)
{
    int i, j;
    double (*d)[src_h] = dest, (*s)[src_w] = src;
    for (i = 0; i < src_h; i++)
        for (j = 0; j < src_w; j++)
            d[j][i] = s[i][j];
}


void print_matrix(int n, int m, double *matrix)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0;j<m;j++)
        {
            printf("%2.15f   ", *((matrix+i*n) + j));
        }
        printf("\n");
    }
}