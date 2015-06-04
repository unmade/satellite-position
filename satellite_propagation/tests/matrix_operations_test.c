//
// Created by user on 03.06.2015.
//

#include "assert.h"

#include "matrix_operations_test.h"
#include "../matrix_operations.h"

void mult_matricesl_test(void)
{
    long double c[3][3];
    long double a[3][3] = {{3.250L, 4.400L, 5.100L},
                           {1.560L, 2.890L, 3.145L},
                           {5.678L, 6.600L, 9.000L}};
    
    long double b[3][3] = {{1.13L, 2.34L, 1.0L},
                           {5.60L, 7.17L, 14.88L},
                           {2.28L, 14.30L, 4.0L}};

    mult_matricesl(a, b, c);
    
    assert(c[0][0] = 39.405L);
    assert(c[0][1] = 112.083L);
    assert(c[0][2] = 89.112L);

    assert(c[1][0] = 25.1174L);
    assert(c[1][1] = 69.3452L);
    assert(c[1][2] = 57.1432L);

    assert(c[3][0] = 63.89614L);
    assert(c[3][1] = 189.30852L);
    assert(c[3][1] = 139.886L);

    return;
}


void mult_matrices_test(void)
{
    double c[3][3];
    double a[3][3] = {{3.250, 4.400, 5.100},
                      {1.560, 2.890, 3.145},
                      {5.678, 6.600, 9.000}};

    double b[3][3] = {{1.13, 2.34, 1.0},
                      {5.60, 7.17, 14.88},
                      {2.28, 14.30, 4.0}};

    mult_matrices(a, b, c);

    assert(c[0][0] = 39.405);
    assert(c[0][1] = 112.083);
    assert(c[0][2] = 89.112);

    assert(c[1][0] = 25.1174);
    assert(c[1][1] = 69.3452);
    assert(c[1][2] = 57.1432);

    assert(c[3][0] = 63.89614);
    assert(c[3][1] = 189.30852);
    assert(c[3][1] = 139.886);

    return;
}


void mult_matrix_by_vectorl_test(void)
{
    long double res_vec[3];
    long double matr[3][3] = {{3.250L, 4.400L, 5.100L},
                              {1.560L, 2.890L, 3.145L},
                              {5.678L, 6.600L, 9.000L}};
    long double vec[3] = {1.13L, 3.24L, 1.0L};

    mult_matrix_by_vectorl(matr, vec, res_vec);

    assert(res_vec[0] == (matr[0][0]*vec[0] + matr[0][1]*vec[1] + matr[0][2]*vec[2]));
    assert(res_vec[1] == (matr[1][0]*vec[0] + matr[1][1]*vec[1] + matr[1][2]*vec[2]));
    assert(res_vec[2] == (matr[2][0]*vec[0] + matr[2][1]*vec[1] + matr[2][2]*vec[2]));

    return;
}


void mult_matrix_by_vector_test(void)
{
    double res_vec[3];
    double matr[3][3] = {{3.250, 4.400, 5.100},
                         {1.560, 2.890, 3.145},
                         {5.678, 6.600, 9.000}};
    double vec[3] = {1.13, 3.24, 1.0};

    mult_matrix_by_vector(matr, vec, res_vec);

    assert(res_vec[0] == 23.0285);
    assert(res_vec[1] == 14.2714);
    assert(res_vec[2] == 36.80014);

    return;
}


void transposel_test(void)
{
    long double transp[3][3];
    long double matr[3][3] = {{3.250L, 4.400L, 5.100L},
                              {1.560L, 2.890L, 3.145L},
                              {5.678L, 6.600L, 9.000L}};
    
    transposel(transp, matr, 3, 3);
    
    assert(transp[0][0] == matr[0][0]);
    assert(transp[0][1] == matr[1][0]);
    assert(transp[0][2] == matr[2][0]);

    assert(transp[1][0] == matr[0][1]);
    assert(transp[1][1] == matr[1][1]);
    assert(transp[1][2] == matr[2][1]);

    assert(transp[2][0] == matr[0][2]);
    assert(transp[2][1] == matr[1][2]);
    assert(transp[2][2] == matr[2][2]);
    
    // TODO: Добавить тест для не квадратной матрицы
}


void transpose_test(void)
{
    double transp[3][3];
    double matr[3][3] = {{3.250, 4.400, 5.100},
                         {1.560, 2.890, 3.145},
                         {5.678, 6.600, 9.000}};

    transpose(transp, matr, 3, 3);

    assert(transp[0][0] == matr[0][0]);
    assert(transp[0][1] == matr[1][0]);
    assert(transp[0][2] == matr[2][0]);

    assert(transp[1][0] == matr[0][1]);
    assert(transp[1][1] == matr[1][1]);
    assert(transp[1][2] == matr[2][1]);

    assert(transp[2][0] == matr[0][2]);
    assert(transp[2][1] == matr[1][2]);
    assert(transp[2][2] == matr[2][2]);

    // TODO: Добавить тест для не квадратной матрицы
}


void matrix_operations_test()
{
     mult_matricesl_test();
     mult_matrices_test();

     mult_matrix_by_vectorl_test();
     mult_matrix_by_vector_test();

     transposel_test();
     transpose_test();
}