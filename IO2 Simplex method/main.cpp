#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <locale>
using namespace std;
typedef double t_base;

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

bool isBV(t_base** a, size_t size1, size_t coln)
{
    size_t i = 0;
    bool f = true, wasOne = false;
    while (i < size1 && f)
    {
        if (a[i][coln] != 0)
        {
            if (a[i][coln] == 1 && !wasOne)
            {
                wasOne = false;
            }
            else
                f = false;
        }
        ++i;
    }
    return f;
}

void printSimplexMat(t_base** a, size_t size1, size_t size2)
{
    size_t sinCol = 0;
    cout << "\nB.v\tFree v\t";
    for (size_t j = 1; j < size2; ++j)
        cout << "x" << j << "\t";
    cout << "\n";
    for (size_t i = 0; i < size1; i++)
    {
        //ñëåäóþùé öèêë ñàì îïðåäåëÿåò êàêàÿ ïåðåìåííàÿ áàçèñíàÿ
        for (size_t j = sinCol+1; j < size2; ++j)
            if (isBV(a, size1, j))
            {
                sinCol = j;
                break;
            }

        if (i == size1-1)
            cout << "z\t";
        else
            cout << "x" << sinCol << "\t";
        //end of self-definition
        for (size_t j = 0; j < size2; j++)
            cout << a[i][j] << '\t';
        cout << endl;
    }
}

void copyMat(t_base** a, t_base** data, size_t sz1, size_t sz2)
{
    for (size_t i = 0; i < sz1; ++i)
        for (size_t j = 0; j < sz2; j++)
            a[i][j] = data[i][j];
}

void swap_str(t_base* a, t_base* b, size_t size)
//îáìåí ìåñòàìè ñòðîê a è b, size - êîë-âî ýëåìåíòîâ â ñòðîêàõ
{
    t_base t;
    for (size_t i = 0; i < size; ++i)
    {
        t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}

bool find_negEl0(t_base* a, size_t sz, size_t& getCol)
{

    size_t j = 0;
    for (; j < sz && a[j] >= 0; ++j);
    if (j < sz)
    {
        getCol = j;
        return 1;
    }
    else
        return 0;
}

bool find_negEl(t_base* a, size_t sz, size_t& getCol)
{

    size_t j = 0;
    bool isNeg = false;
    while (j < sz && !isNeg)
    {
        isNeg = a[j] < 0;
        if (isNeg)
        {
            size_t i_maxNeg = j++;
            while (j < sz)
            {
                if (a[j] < 0 && fabs(a[j]) > fabs(a[i_maxNeg]))
                    i_maxNeg = j;
                ++j;
            }
            getCol = i_maxNeg;
            break;
        }
        ++j;
    }
    if (isNeg)
        return 1;
    else
        return 0;
}

bool findPositCol(t_base** a, size_t colNum, size_t sz1, size_t& getRow)
{
    size_t i = 0;
    for (; i < sz1 && a[i][colNum] < 0; ++i);
    if (i < sz1)
    {
        getRow = i;
        return 1;
    }
    else
        return 0;
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
        printSimplexMat(a, size1, size2);
    }
}

void SimplexMethod(t_base** a, size_t sz1, size_t sz2)
{
    size_t leadRow, leadCol;
    while (find_negEl(a[sz1-1] + 1, sz2 - 1, leadCol))
    {
        leadCol++;
        //îòáîð âåäóùåãî ýëåìåíòà
        if (findPositCol(a, leadCol, sz1, leadRow))
        {
            t_base min = a[leadRow][leadCol] / a[leadRow][0];
            size_t i_min = leadRow;
            while (i_min < sz1)
            {
                if (a[i_min][leadCol] >= 0
                     && a[i_min][leadCol] / a[i_min][0] < min)
                {
                    min = a[i_min][leadCol] / a[i_min][0];
                    leadRow = i_min;
                }
                i_min++;
            }
        //êîíåö îòáîðà âåäóùåãî ýëåìåíòà
            makeSingularColumn(a, sz1, sz2, leadRow, leadCol);
        }
        else
        {
            cout << "There is no solutions\n";
            return;
        }
    }
}

void fMatInput(char* f_name, t_base** a, size_t size1, size_t size2)
//ââîä èç ôàéëà ñ èì f_name ìàòð à ðàçìåðà size2 * size2
{
    ifstream f(f_name);
    for (size_t i = 0; i < size1; ++i)
        for (size_t j = 0; j < size2; ++j)
            f >> a[i][j];
    f.close();
}

void printSolution(t_base** a, size_t sz1, size_t sz2)
{
    size_t sinCol = 0;
    cout << "zmax = fmin = " << a[sz1-1][0] << "\n";
    for (size_t i = 0; i < sz1-1; i++)
    {
        for (size_t j = sinCol+1; j < sz2; ++j)
                if (isBV(a, sz1, j))
                {
                    sinCol = j;
                    break;
                }
        cout << "x" << sinCol << " = " << a[i][0] << "\t";

    }
}

int main()
{
    setlocale(LC_ALL, "Rus");
    size_t sz1, sz2;
    char ckey;

    cout << "input the sizes that cover whole simplex table data\n";
    cin >> sz1 >> sz2;

    cout << "Do you want to input data via console or file\n\
    press f to input via file or another key to input via console\n";
    cin >> ckey;

    t_base** a = newMat(sz1, sz2);

    if (ckey == 'f')
    {
        cout << "input the simplex matrix in sxslae.txt\n";
        char f_name[20] = "sxslae.txt";
        fMatInput(f_name, a, sz1, sz2);
    }
    else
    {
        cout << "Введите симплекс-таблицу\n";
        for (size_t i = 0; i < sz1; ++i)
            for (size_t j = 0; j < sz2; ++j)
                cin >> a[i][j];
    }

   // cout << "1st simplex table";

    printSimplexMat(a, sz1, sz2);

    SimplexMethod(a, sz1, sz2);
   // cout << "Last Simplex-table";
    //printSimplexMat(a, sz1, sz2);
    //cout << "Результат\n";
    printSolution(a, sz1, sz2);
    system("pause");
    return 0;
}
