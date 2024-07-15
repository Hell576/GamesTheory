#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"


t_base** alloc(unsigned n, unsigned m)
//возвр выделенную в памяти матрицу m * n
{
    t_base **pa = (t_base**)calloc(n, sizeof(t_base*));
    for (unsigned i = 0; i < n; ++i)
    {
        pa[i] = (t_base*)calloc(m, sizeof(t_base));
    }
    return pa;
}

void output(t_base** a, size_t size1, size_t size2)
//вывод матр а размера size2 * size2
{
    printf("\n");
    for (size_t i = 0; i < size1; ++i)
    {
        for (size_t j = 0; j < size2; ++j)
            printf("%f ", a[i][j]);
        printf("\n");
    }
}

void fmat_input(char* f_name, t_base** a, size_t size1, size_t size2)
//ввод из файла с им f_name матр а размера size2 * size2
{
    FILE* f = fopen(f_name, "r");
    for (size_t i = 0; i < size1; ++i)
        for (size_t j = 0; j < size2; ++j)
            fscanf(f, "%f", &a[i][j]);
    fclose(f);
}

void f_slae_input(char* f_name, t_base** a, size_t size1,
                  t_base** p, size_t colp)
//ввод слау. а - матрица коэфф при неизвестных размера size1*size2
//p - матрица свободных коэфф размера size1 * colp
{
    FILE* f = fopen(f_name, "r");
    for (size_t i = 0; i < size1; ++i)
    {
        for (size_t j = 0; j < size1; ++j)
            fscanf(f, "%f", &a[i][j]);
        for (size_t j = 0; j < colp; ++j)
            fscanf(f, "%f", &p[i][j]);
    }
    fclose(f);
}
void trans_mat(t_base** a, size_t n)
//транспонирование матрицы a размера n * n
{
    t_base t;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            t = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = t;
        }
    }
}

void mult_mat(const t_base** a, const t_base** b, t_base** c,
              size_t n1, size_t m2, size_t k)
//умножение матрицы а на матрицу b, n1 - кол-во строк а,
//m2 кол-во столбов b, k - кол-во элементов в строках а и столбцах b
{
    for(size_t i = 0; i < n1; ++i)
    {
        for (size_t j = 0; j < m2; ++j)
        {
            int res = 0;
            for (size_t ij = 0; ij < k; ++ij)
            {
                res += a[i][ij] * b[ij][j];
            }
            c[i][j] = res;
        }
    }
}

void n_m_mat(const t_base** a, t_base** b, size_t size1, size_t size2, t_base k)
//умножение матрицы а размера size1 * size2 на число k, b - результат умножения матрицы
{
    for (size_t i = 0; i < size1; ++i)
    {
        for (size_t j = 0; j < size2; ++j)
        {
            b[i][j] = a[i][j] * k;
        }
    }
}

void plus_mat(const t_base** a, const t_base** b, t_base** c, size_t n, size_t m)
//сложение матриц а и b, обе размера n * m, c - результат сложения матриц
{
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
}

void minus_mat(const t_base** a, const t_base** b, t_base** c, size_t n, size_t m)
//вычитание матриц а и b, обе размера n * m, c - результат вычитания матриц
{
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
}

void n_mult_str(t_base n, const t_base* a, t_base* b, size_t size)
//умножение строки(числового вектора) а размера size на число n, b - результат умножения
{
    for (size_t i = 0; i < size; ++i)
    {
        b[i] = a[i] * n;
    }
}

void plus_str(t_base* a, const t_base* b, size_t size)
//прибавление к вектору а размера size вектор b
{
    for (size_t i = 0; i < size; ++i)
    {
        a[i] = a[i] + b[i];
    }
}

t_base mmax(t_base** a, size_t size, unsigned j, unsigned* m_str)
//поиск максимального по модулю элемента в столбце j матрицы а,
// где size - кол-во строк в матрице, m_str - номер строки с максимальным элементом
{
    t_base max = (*(a + j)) [j];
    *m_str = j;
    for (size_t i = j + 1; i < size; ++i)
    {
        if ( fabs( (*(a + i))[j] ) > fabs(max) )
        {
            max = (*(a + i))[j];
            *m_str = i;
        }
    }
    return max;
}

void swap_var(t_base* a, t_base* b)
{
    t_base t = *a;
    *a = -(*b);
    *b = -t;
}

void swap_str(t_base* a, t_base* b, size_t size)
//обмен местами строк a и b, size - кол-во элементов в строках
{
    t_base t;
    for (size_t i = 0; i < size; ++i)
    {
        t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}


t_base sum(t_base* a, t_base* x, size_t size, size_t k)
//сложение элементов строки(вектора) а размера size, взятой(го) из матрицы,
//x - вектор неизвестных, k - номер строки матрицы.
{
    t_base res = 0;
        for (size_t j = k + 1; j < size; j++)
        {
            res += a[j] * x[j];
        }
    return res;
}

int forward_stroke(t_base** a, size_t size1, size_t size2, t_base** p, size_t colp)
//прямой ход метода Гаусса. a - матрица коэфф при неизвестных размера size1*size2,
//colp - кол-во столбцов матрицы свободных коэфф p, x - вектор неизвестных

{
   // t_base multRow[size2 + colp];
    t_base* multRow = (t_base*)calloc(size2 + colp, sizeof(t_base));
    int k = 0;
        for (size_t j = 0; j < size2 - 1; j++)
        {
            size_t lead_str;
            t_base leadEl = -1 * mmax(a, size1, j, &lead_str);

                if ( (j != lead_str) || (a[j][0] != 1) )
                {
                    swap_str(a[j], a[lead_str], size2);
                //j can go either through the lines and column
                    k++;
                        if (p != NULL)
                            swap_str(p[j], p[lead_str], colp);
                }

                if (leadEl != 0)
                {
                    for (size_t i = j + 1; i < size1; ++i)
                    {
                    t_base mu = a[i][j] / leadEl;
                    n_mult_str(mu, a[j], multRow, size2);
                    plus_str(a[i], multRow, size2);
                    if (p != NULL)
                    {
                        n_mult_str(mu, p[j], multRow + size2, colp);
                        plus_str(p[i], multRow + size2, colp);
                    }
                    }
                }
        }
    free(multRow);
    return k;
}

t_base det(t_base** a, size_t size1)
//определитель матрицы а размера size1*size1
{
    int k = forward_stroke(a, size1, size1, NULL, 0);
    t_base det = 1;
        for (size_t i = 0; i < size1; i++)
        {
            det *= a[i][i];
        }
    return pow(-1, k) * det;
}

int backward_stroke(t_base** a, t_base** p1, int j, t_base* x,
                      size_t size1, size_t size2)
//обратный ход метода Гаусса. a - матрица коэфф при неизвестных размера size1*size2,
//j - столбец матрицы свободных коэфф p1, x - вектор неизвестных

{
    int flag = 0;
        if (a[size1-1][size2-1] != 0)
            x[size1 -1] = p1[size1 -1][j] / a[size1 -1][size2 -1];
        else if (p1[size1 -1] == 0)
                flag = 1;
                //x[size1 -1] = 0;
             else return -1;

        for (size_t k = size1 - 2; k >= 0; k--)
        {
                if (flag != 1)
                {

                }
            t_base s = sum(a[k], x, size2, k);
                if (s != 0)
                    x[k] = (p1[k][j] - s) / a[k][k];
                else if (p1[k][j] == 0)
                        x[k] = 0;
                     else return -1;
        }
    return 0;
}

void release_mem(t_base **pa, size_t n)
//освобождение матрицы по адрессу ра, n - кол-во строк матрицы
{
    for (size_t i = 0; i < n ; ++i)
    {
      free(pa[i]);
    }
  free(pa);
}
