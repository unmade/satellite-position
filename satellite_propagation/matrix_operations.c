#include <stdio.h>
#include "matrix_operations.h"


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


//void print_matrix( int n, int m, double **matrix)
//{
//    for (int i=0; i<n; i++)
//    {
//        for (int j=0;j<m;j++)
//        {
//            printf("%2.15f   ", matrix[i][j]);
//        }
//        printf("\n");
//    }
//}