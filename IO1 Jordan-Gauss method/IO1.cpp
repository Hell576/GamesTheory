#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "matrix.h"
using namespace std;

const float eps = 0.00001;

bool isZero(t_base x)
{
    return fabs(x) < eps;
}

bool isZeroRow(t_base* a, size_t sz)
{
    size_t i;
    for (i = 0; i < sz && isZero(a[i]); ++i);
    return i == sz;
}

t_base** newMat(size_t sz1, size_t sz2)
{
    t_base** a = new t_base*[sz1];
    for (size_t i = 0; i < sz1; ++i)
        a[i] = new t_base[sz2];
    return a;
}

void delMat(t_base** a, size_t sz1)
{
    for (size_t i = 0; i < sz1; ++i)
        delete[] a[i];
    delete[] a;
}

void printMat(t_base** a, size_t size1, size_t size2)
{
    for (size_t i = 0; i < size1; i++)
    {
        if (!isZeroRow(a[i], size2))
        {
            for (size_t j = 0; j < size2; j++)
                cout << a[i][j] << '\t';
            cout << endl;
        }
    }
    cout << endl;
}

void printExtMat(t_base** a, size_t size1, size_t size2)
{
    for (size_t i = 0; i < size1; i++)
    {
        if (!isZeroRow(a[i], size2))
        {
            for (size_t j = 0; j < size2; j++)
            {
                if (a[i][j] == 0)
                    a[i][j] = 0;
                cout << a[i][j];
                if (j == size2-2)
                    cout << '|';
                cout << "\t";
            }
            cout << endl;
        }
    }
    cout << endl;
}

void copyMat(t_base** a, t_base** data, size_t sz1, size_t sz2)
{
    for (size_t i = 0; i < sz1; ++i)
        for (size_t j = 0; j < sz2; j++)
            a[i][j] = data[i][j];
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

bool find_nonzero(t_base** a, size_t numRow, size_t col, size_t& row)
//поиск ненулевого элемента в столбце col матрицы а, возврат 1, если есть, иначе - 0.
//где numRow - кол-во строк в матрице, a row - номер строки ненулевого элемента
{
    for (row = 0; row < numRow && ((*(a + row))[col] == 0); ++row);
    if (row < numRow)
        return true;
    else return false;
}
void makeSingularColumn(t_base** a, size_t size1, size_t size2,
                        size_t iLead, size_t jLead)
{
    t_base leadEl = a[iLead][jLead];
    if (!isZero(leadEl))
    {
        for (unsigned j = 0; j < size2; j++)
            a[iLead][j] /= leadEl;
        for (unsigned i = 0; i < size1; ++i)
        {
            if (i != iLead)
            {
                t_base k = a[i][jLead];
                for (unsigned j = 0; j < size2; j++)
                {
                    a[i][j] += -1 * k * a[iLead][j];
                }
            }
        }
        //printExtMat(a, size1, size2);
    }
}

size_t Jordan_Gauss(t_base** a, unsigned* basisVars, size_t size1, size_t size2)
{
    size_t suitRow, rank = size1;
    for (unsigned i = 0; i < size1; ++i)
    {

        if (!isZero(a[i][basisVars[i] - 1]))
        {
            makeSingularColumn(a, size1, size2, i, basisVars[i]-1);
        }
        else //если ведущий елемент = 0
        {
            if (isZeroRow(a[i], size2)) //если строка нулевая или
                rank--;
            else if (isZeroRow(a[i], size2-1) && a[size2-1] != 0)///!!!!
            {
                cerr << "The matrix has no solves\n";
                system("pause");
                exit(-1);
            }
            //или только ведущий элем = 0, мы будем искать строку с ненулевым вед элемент
            //и при возможности обменяем негодную строку на более подходящую
            if (find_nonzero(a + i, size1 - i, i, suitRow))
                swap_str(a[i], a[suitRow], size2);
        }
    }
    return rank;
}

bool isOpornBasisSol(t_base** mat, size_t sz1, size_t sz2)
{
	bool ans = 1;
	for (size_t i = 0; i < sz1 && ans; i++)
		ans = mat[i][sz2-1] >= eps;
	return ans;
}

typedef t_base aimFun(t_base*, size_t);

t_base z_fun(t_base* x, size_t sz)
{
    return 3 * x[0] + 0.5 * x[1];
}

void max_f(t_base* a, t_base* maxVars, size_t sz, aimFun f)
{
    if (f(a, sz) > f(maxVars, sz))
        for (size_t i = 0; i < sz; ++i)
            maxVars[i] = a[i];
}

//Вывод опорного базисного решения на экран
void printOpSol(t_base** mat, unsigned* arrans, size_t sz1, size_t sz2)
{
	for (size_t i = 0; i < sz1; i++)
		cout << "x" << arrans[i] << " is " << mat[i][sz2-1] << "\n";
}

size_t iterateBasisVars(unsigned* basisVars, t_base* maxVars, unsigned minStartEl, unsigned i, unsigned n, unsigned k,
                  t_base** a, t_base **temp)
{
    //don't understand what is minStartEl? put 1 as a parameter

    static size_t q;
    for (unsigned x = minStartEl; x < n - k + i + 2; x++)///it works don't touch it!
    {
        unsigned col = n+1;
        basisVars[i] = x;
        if (i == k - 1)
        {
            for (size_t o = 0; o < k; o++)
                cout << "x" << basisVars[o] << ' ';
            cout << endl;

            copyMat(temp, a, k, col);
            Jordan_Gauss(temp, basisVars, k, col);
            printExtMat(temp, k, col);

            q++;

            if (isOpornBasisSol(temp, k, col))
            {
                cout << "reference solution\n";
                t_base currVars[k];
                for (size_t o = 0; o < k; o++)
                    currVars[o] = temp[o][n];

                printOpSol(temp, basisVars, k, col);
                max_f(currVars, maxVars, k, z_fun);

            }
        }
        else
        {
            iterateBasisVars(basisVars, maxVars, x + 1, i + 1, n, k,
                       a, temp);
        }
    }
    return q;
}

void Sochet(t_base** a, unsigned n, unsigned k)//hardcode
{
    //don't understand what is minStartEl? put 1 as a parameter
    unsigned* c = new unsigned[k]; t_base* maxVars = new t_base[k]{};
    t_base **temp = newMat(k, n+1);

    iterateBasisVars(c, maxVars, 1, 0, n, k, a, temp);

    //delMat(temp, k);
    cout << "max fun val is " << z_fun(maxVars, k) << endl;
    cout << "the vars are " << maxVars[0] << " " << maxVars[1] << "\n";
    delete[] c; delete[] maxVars;
}

void fMatInput(char* f_name, t_base** a, size_t size1, size_t size2)
//ввод из файла с им f_name матр а размера size2 * size2
{
    ifstream f(f_name);
    for (size_t i = 0; i < size1; ++i)
        for (size_t j = 0; j < size2; ++j)
            f >> a[i][j];
    f.close();
}



/*bool isBasisSolves(t_base** a, size_t sz1, size_t sz2)
{
    size_t i, j = sz2 - 1;
    for (i = 0; i < sz1 && a[i][j] == 0; i++);
    if (i < sz1)
        return false;
    else return true;
} */



int main()
{
    unsigned n, m;

 /*   cout << "the processed matrix must be square & will be in the project dir\nin file \"slae.txt\" \
,follow the format:\na11 a12 ... a1N b1\n\
a21 a22 ... a2N b2\n.\n.\n.\naM1 aM2 ... aMN bM1 bM2 ... bN\n";
 */
 //   cout << "Are you ready? ";
 cout << "z = 3*x1 + 0.5*x2\n";
   // system("pause");

    char f_name[] = "slae.txt";
    cout << "Input #n rows & #m columns(including coloums for b)\n";
    cin >> n >> m;
    unsigned vars = m - 1;
    t_base** a = newMat(n, m);

    fMatInput(f_name, a, n, m);
    cout << "source matrix\n";
    printExtMat(a, n, m);
    //решаем слау методом Жордана-Гаусса обычным способом
    unsigned bv[vars];
    for (unsigned i = 0; i < vars; ++i)
        bv[i] = i+1;

    t_base **temp = newMat(n, m);
    copyMat(temp, a, n, m);

    size_t rank = Jordan_Gauss(temp, bv, n, m);
    delMat(temp, n);
    //ищем множество базисных видов(решений) перебором
    cout.precision(2);//выводить 2 цифры после запятой
    Sochet(a, vars, rank);

    delMat(a, n);

    system("pause");
    return 0;
}
