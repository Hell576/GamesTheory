#ifndef MATRIX_H
#define MATRIX_H
typedef double t_base;

t_base** alloc(unsigned n, unsigned m);
//����� ���������� � ������ ������� m * n



void output(t_base** a, size_t size1, size_t size2);
//����� ���� � ������� size2 * size2



void fmat_input(char* f_name, t_base** a, size_t size1, size_t size2);
//���� �� ����� � �� f_name ���� � ������� size2 * size2



void f_slae_input(char* f_name, t_base** a, size_t size1,
                   t_base** p, size_t colp);
//���� ����. � - ������� ����� ��� ����������� ������� size1*size2
//p - ������� ��������� ����� ������� size1 * colp



void trans_mat(t_base** a, size_t n);
//���������������� ������� a ������� n * n



void mult_mat(const t_base** a, const t_base** b, t_base** c,
              size_t n1, size_t m2, size_t k);
//��������� ������� � �� ������� b, n1 - ���-�� ����� �,
//m2 ���-�� ������� b, k - ���-�� ��������� � ������� � � �������� b



void n_m_mat(const t_base** a, t_base** b,
             size_t size1, size_t size2, t_base k);
//��������� ������� � ������� size1 * size2 �� ����� k, b - ��������� ��������� �������



void n_mult_str(t_base n, const t_base* a, t_base* b, size_t size);
//��������� ������(��������� �������) � ������� size �� ����� n, b - ��������� ���������



void plus_str(t_base* a, const t_base* b, size_t size);
//����������� � ������� � ������� size ������ b



t_base mmax(t_base** a, size_t size, unsigned j, unsigned* m_str);
//����� ������������� �� ������ �������� � ������� j ������� �,
// ��� size - ���-�� ����� � �������, m_str - ����� ������ � ������������ ���������



void swap_var(t_base* a, t_base* b);
//����� ������� ���������� a � b



void swap_str(t_base* a, t_base* b, size_t size);
//����� ������� ����� a � b, size - ���-�� ��������� � �������



void release_mem(t_base **pa, size_t n);
//������������ ������� �� ������� ��, n - ���-�� ����� �������

#endif // MATRIX_H
