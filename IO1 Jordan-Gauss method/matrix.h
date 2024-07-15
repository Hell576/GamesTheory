#ifndef MATRIX_H
#define MATRIX_H
typedef double t_base;

t_base** alloc(unsigned n, unsigned m);
//возвр выделенную в памяти матрицу m * n



void output(t_base** a, size_t size1, size_t size2);
//вывод матр а размера size2 * size2



void fmat_input(char* f_name, t_base** a, size_t size1, size_t size2);
//ввод из файла с им f_name матр а размера size2 * size2



void f_slae_input(char* f_name, t_base** a, size_t size1,
                   t_base** p, size_t colp);
//ввод слау. а - матрица коэфф при неизвестных размера size1*size2
//p - матрица свободных коэфф размера size1 * colp



void trans_mat(t_base** a, size_t n);
//транспонирование матрицы a размера n * n



void mult_mat(const t_base** a, const t_base** b, t_base** c,
              size_t n1, size_t m2, size_t k);
//умножение матрицы а на матрицу b, n1 - кол-во строк а,
//m2 кол-во столбов b, k - кол-во элементов в строках а и столбцах b



void n_m_mat(const t_base** a, t_base** b,
             size_t size1, size_t size2, t_base k);
//умножение матрицы а размера size1 * size2 на число k, b - результат умножения матрицы



void n_mult_str(t_base n, const t_base* a, t_base* b, size_t size);
//умножение строки(числового вектора) а размера size на число n, b - результат умножения



void plus_str(t_base* a, const t_base* b, size_t size);
//прибавление к вектору а размера size вектор b



t_base mmax(t_base** a, size_t size, unsigned j, unsigned* m_str);
//поиск максимального по модулю элемента в столбце j матрицы а,
// где size - кол-во строк в матрице, m_str - номер строки с максимальным элементом



void swap_var(t_base* a, t_base* b);
//обмен местами переменных a и b



void swap_str(t_base* a, t_base* b, size_t size);
//обмен местами строк a и b, size - кол-во элементов в строках



void release_mem(t_base **pa, size_t n);
//освобождение матрицы по адрессу ра, n - кол-во строк матрицы

#endif // MATRIX_H
