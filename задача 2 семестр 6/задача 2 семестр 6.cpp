// задача 2 семестр 6.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
using namespace std;
//коэффициенты диф.ур.
double P(double x)
{
    return (1) / (2 + (x / 3));
}
double Q(double x)
{
    return (exp(x / 5));
}
//Копирует матрицу B в матрицу A
void cop_mat(double** A, double** B, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = B[i][j];
}
//Копирует вектор b в вектор a
void cop_vect(double* a, double* b, int n)
{
    for (int i = 0; i < n; i++)
        a[i] = b[i];
}
//выводит матриицу A порядка n
void show_mat(double** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << fixed << setprecision(6) << " \t" << A[i][j] << " \t";
        }
        std::cout << endl;
    }
}
//выводит вектор V  
void show_vect(double* V, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << fixed << setprecision(6) << " \t" << V[i] << " \t";
        std::cout << endl;
    }
}
//перемножение матриц
double** mul_mat(double** A, double** B, int n)
{
    double** C = new double* [n];
    for (int i = 0; i < n; i++)
    {
        C[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
            //cout << C[i][j] << " ";
        }
        //cout << endl;
    }
    // show(C,n);
    return C;
}
//Умножение матрицы на вектор
double* Mat_x_vect(double** Mat, double* vect, int n)
{
    double* res = new double[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
        for (int j = 0; j < n; j++)
        {
            res[i] += Mat[i][j] * vect[j];
        }
    }
    return res;
}
//математическая функция sign(x)
double sign(double value)
{
    if (value > 0)
    {
        return 1;
    }
    else if (value < 0)
    {
        return -1;
    }
    else if (value == 0)
    {
        return 0;
    }
}
//строит нулевой фектор порядка n
double* get_null_vect(int n)
{
    double* V = new double[2 * n];
    for (int i = 0; i <= n; i++)
    {
        V[i] = 0.0;
    }
    return V;
}
//строит нулевую матрицу порядка n 
double** get_null_mat(int n)
{
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}
//транспонирует матрицу А порядка n
double** transp_mat(double** A, int n)
{
    double tmp;
    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
        {
            tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    return A;
}
//для решения системы с матрицей "a" вектором правых частей "y" и порядком "n"
double* gauss2(double** a, double* y, int n)
{
    double* x, max;
    int k, index;
    const double eps = 0.00001;  // точность
    x = new double[n];
    k = 0;
    for (k = 0; k < n; k++)
    {
        // Поиск строки с максимальным a[i][k]
        max = fabs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(a[i][k]) > max)
            {
                max = fabs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (fabs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}
//Многочлены Лежандра
double Legendre(int n, double x)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return x;
    else
        return (Legendre(n - 1, x) * x * (2 * n - 1) / (n)) - (Legendre(n - 2, x) * (n - 1) / (n));
}
//Локализация корней многочленов Лежандра
double* loc_root(double a, double b, double N, int n)
{
    double* mas = new double[100];
    mas[0] = n;
    double h = (b - a) / N;
    double tmp = 0;
    int j = 2;
    while (a + h <= b)
    {
        if ((Legendre(n, a) * Legendre(n, a + h)) <= 0)
        {
            mas[j] = a;
            mas[j + 1] = a + h;
            j += 2;
        }
        a = a + h;
    }
    if (mas[0] > 0)
    {
        //cout << fixed << setprecision(0) << "количество корней ортоганального многочлена на отрезке (a,b) = " << mas[0] << endl;
        return mas;
    }
    else
    {
        //cout << "нет корней на данном отрезке" << endl;
    }
    return 0;
}
//Коэффициенты для формулы Гаусса
double* koeff(double* z, int n)
{
    double* A = new double[100];
    for (int i = 0; i < n; i++)
    {
        A[i] = (2 * (1 - pow(z[i], 2)) / (pow(n, 2) * Legendre(n - 1, z[i]) * Legendre(n - 1, z[i])));
    }
    return A;
}
//Поиск корней многочленов Лежандра
double Secant(double a, double b, double EPS, int n)
{
    double fapprox, sapprox, x;
    int counter = 0;

    fapprox = a;
    sapprox = b;

    x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
    counter++;

    while (abs(x - sapprox) >= EPS)
    {
        fapprox = sapprox;
        sapprox = x;
        x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
        counter++;
    }

    return x;
}
//находит максимальный верхнедиагональный элемент матрицы А
double* Max_upper_diag(double** A, int n)
{
    double* res = get_null_vect(n);
    res[0] = 0.0;
    res[1] = 0.0;
    res[2] = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                if (i < j)
                {
                    if (fabs(A[i][j]) > res[0])
                    {
                        res[0] = fabs(A[i][j]);
                        res[1] = i;
                        res[2] = j;
                    }
                }
            }
        }
    }
    return res;
}
//готовит две матрицы и вектор для шага метода вращений 
void Prepare(double** X, double** A_next, double** V, int n)
{
    for (int i = 0; i < n; i++)
    {
        X[i] = new double[n];
        A_next[i] = new double[n];
        V[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                X[i][j] = 1;
                A_next[i][j] = 0.0;
                V[i][j] = 1;
            }
            else
            {
                X[i][j] = 0.0;
                A_next[i][j] = 0.0;
                V[i][j] = 0.0;
            }
        }
    }
}
//Метод вращений Якоби
double* Jacob(double** A, int n)
{
    double EPS = 0.000001;
    double** X = new double* [n];
    double** A_next = new double* [n];
    double** V = new double* [n];
    double c = 0.0;
    double s = 0.0;
    Prepare(X, A_next, V, n);
    double* MAS = Max_upper_diag(A, n);
    double max = MAS[0];
    int i_k = MAS[1];
    int j_k = MAS[2];
    while (max >= EPS)
    {
        double d = sqrt(pow((A[i_k][i_k] - A[j_k][j_k]), 2) + 4 * pow(A[i_k][j_k], 2));
        c = sqrt(0.5 * (fabs((A[i_k][i_k] - A[j_k][j_k])) / d) + 0.5);
        s = sign(A[i_k][j_k] * (A[i_k][i_k] - A[j_k][j_k])) * (sqrt(-0.5 * (fabs((A[i_k][i_k] - A[j_k][j_k])) / d) + 0.5));
        V[i_k][i_k] = c;
        V[j_k][j_k] = c;
        V[i_k][j_k] = (-1) * s;
        V[j_k][i_k] = s;
        cop_mat(X, mul_mat(X, V, n), n);
        A_next = mul_mat(transp_mat(V, n), mul_mat(A, V, n), n);
        cop_mat(A, A_next, n);
        Prepare(V, A_next, V, n);
        double* MAS_next = Max_upper_diag(A, n);
        max = MAS_next[0];
        i_k = MAS_next[1];
        j_k = MAS_next[2];
    }
    double* eigenvalues = new double[n];
    std::cout << "Собственные значения (Метод Вращения): " << endl;
    for (int i = 0; i < n; i++)
    {
        eigenvalues[i] = A[i][i];
        std::cout << eigenvalues[i] << endl;
    }
    std::cout << "Матрица собственныйх веркторв X: " << endl;
    show_mat(X, n);
    std::cout << endl;
    std::cout << "Проверка длины собственных векторов: " << endl;
    double* length = new double[n];
    for (int i = 0; i < n; i++)
    {
        length[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            length[i] += X[i][j] * X[i][j];
        }
        std::cout << "        " << sqrt(length[i]) << "       ";
    }
    cout << endl;
    return eigenvalues;
}
//строит массив значений системы многочленов Лежандра в данной точке (в данной задаче k=0)
double* Legendre_sys(int k, int n, double x)
{
    double* P_x = get_null_vect(2 * n);
    double a, b, c, d = 0.0;
    P_x[0] = 1;
    P_x[1] = (1 + k) * x;
    for (int i = 0; i < n - 2; i++)
    {
        a = (i + k + 2);
        b = (2 * i + 2 * k + 3);
        c = (i + k + 1);
        d = (i + 2 * k + 2) * (i + 2);
        // cout<<k<<" --- "<<a<<" ---- "<<b<<" ---- " << c << " --- "<<d<<endl;
        P_x[i + 2] = (a * b * x * P_x[i + 1] - a * c * P_x[i]) / d;
    }
/*
        cout << endl;
        cout << "Массив зачений системы многочленов Лежандра n = " << n <<" k =  "<<k<<"  в точке = " << x << endl;
        for (int i = 0; i < n; i++)
        {
            cout <<i<<"   "<< P_x[i] << " "<<endl;
        }
*/  
    return P_x;
}
//строит массив значений координатоной системы многочленов Лежандра в данной точке (узле для вычисление интеграла по Гауссу)
double* Legendre_coord_sys(int n, double x)
{
    double* P_x = Legendre_sys(0,n, x);
    double* phi = get_null_vect(2 * n);
    phi[0] = 1.0 * sqrt(0.5);
    phi[1] = x * sqrt(1.5);
    for (int i = 3; i <= n; i++)
    {
        phi[i - 1] = sqrt(i - 0.5) * P_x[i - 1];
    }
 /*
        cout << endl;
        cout << "Массив зачений координатой системы многочленов Лежандра n = " << n << " в точке = " << x << endl;
        for (int i = 0; i < n; i++)
        {
            cout << i << "   " << phi[i] << " " << endl;
        }
        cout << endl;
*/
    return phi;
}
//строит массив значений производных координатоной системы многочленов лежандра в данной точке (узле для вычисление интеграла по Гауссу)
double* Legendre_Dir_coord_sys(int n, double x)
{
    double* P_x = Legendre_sys(1, n, x);
    double* phi_dir = get_null_vect(2 * n);
    double a, b, c, d = 0.0;
    double e = (1 - pow(x, 2));
    phi_dir[0] = 0.0;
    phi_dir[1] = sqrt(1.5);
    for (int i = 3; i <= n; i++)
    {
        a = sqrt(i - 0.5);
        b = i / 2.0;
        phi_dir[i - 1] = a * b * P_x[i-2];
    }
/*
        cout << endl;
        cout << "Массив зачений производных координатой системы многочленов Лежандра n = " << n << " в точке = " << x << endl;
        for (int i = 0; i < n; i++)
        {
            cout << i << "   " << phi_dir[i] << " " << endl;
        }
*/
    return phi_dir;
}
//строит массив значений вторых производных координатоной системы многочленов лежандра в данной точке (узле для вычисление интеграла по Гауссу)
double* Legendre_Sec_Dir_coord_sys(int n, double x)
{
    double* P_x = Legendre_sys(2, n, x);
    double* phi_sec_dir = get_null_vect(2 * n);
    double a, b, c, d = 0.0;
    double e = (1 - pow(x, 2));
    phi_sec_dir[0] = 0;
    phi_sec_dir[1] = 0;
    for (int i = 3; i <= n; i++)
    {
        a = sqrt(i - 0.5);
        b = i / 2.0;
        c = (i + 1) / 2.0;
//        cout <<i<<"   " << a << "   " << b << "   " << c << "   " << P_x[i - 3] << endl;
        phi_sec_dir[i - 1] = a * b * c * P_x[i - 3];
    }
/*/
    cout << endl;
    cout << "Массив зачений вторых производных координатой системы многочленов Лежандра n = " << n << " в точке = " << x << endl;
    for (int i = 0; i < n; i++)
    {
        cout << i << "   " << phi_sec_dir[i] << " " << endl;
    }
*/
    return phi_sec_dir;
}
//Получает матрицу G_L из метода Ритца
double** Ritz_G_L(int a, int b, int n)
{
    double alpha1 = -0.8;//коэффициенты кравевых условий 
    double alpha2 = -1;
    double beta1 = 0.85;
    double beta2 = 1;
    double aa = alpha1 / alpha2;
    double bb = beta1 / beta2;
    double** phi = get_null_mat(n);
    double** phi_dir = get_null_mat(n);
    double** phi_sec_dir = get_null_mat(n);
    double* z = get_null_vect(2 * n);        //массив для узлов для формулы Гаусса
    int val = 8;                             // количество узлов для формулы Гаусса
    double* MASS = loc_root(-1, 1, 100, val);//локализуем узлы многочлена Лежандра
    double counter = MASS[0];
    int cnt = 0;
    for (int j = 2; j <= 2 * counter + 1; j += 2)
    {
        z[cnt] = Secant(MASS[j], MASS[j + 1], 0.0000000001, val);//находим узлы (корни многочлена Лежандра)
        cnt++;
    }
    if (val % 2 > 0)
    {
        int m = val / 2;
        z[m] = 0.0;
    }
    double* A = koeff(z, val); //находим коэффициенты формулы гаусса
// cout << "Узлы и коэффициенты для формулы Гаусса: " << endl;
// cout << "      Узел      " << " <-> " << "    Коэффициент    " << endl;
    for (int k = 0; k < cnt; k++)
    {
        z[k] = (b - a) / 2 * z[k] + (b + a) / 2; //узлы для формулы гаусса
        A[k] = (b - a) / 2 * A[k];
        // cout << fixed << setprecision(2) << z[k] << "            <->            " << A[k] << endl;
    }
//получим массивы значений координатной системы и ее производной в узлах для формулы гаусса и запишем их в матрицы, чтобы потом из них собирать значения подыитегральной функции.
    double intg = 0.0;
//cout << "Найдем массивы значений функций из координатой системы и их производных в узлах формулы Гаусса: " << endl;
    for (int k = 0; k < cnt; k++)
    {
        phi[k] = Legendre_coord_sys(n, z[k]);
        phi_dir[k] = Legendre_Dir_coord_sys(n, z[k]);
        phi_sec_dir[k] = Legendre_Sec_Dir_coord_sys(n, z[k]);
    }
    phi[cnt] = Legendre_coord_sys(n, a);
    phi[cnt+1] = Legendre_coord_sys(n, b);
    /*
        cout << "Соберем их в матрицы, где i-я строка - значения производной функций на i-ом узле формулы Гаусса:" << endl;
        show(phi, n);
        cout << "------------------------------------------------------------------------------------------------" << endl;
        show(phi_dir, n);
        cout << "------------------------------------------------------------------------------------------------" << endl;
        show(phi_sec_dir, n);
        cout << "------------------------------------------------------------------------------------------------" << endl;

        cout << "Используя эти матрицы найдем коэффициенты исходной системы как интегралы, используя узлы и коэффициенты формулы Гаусса:" << endl;
    */
    double** G = get_null_mat(n);
    double** G_L = get_null_mat(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            
            for (int k = 0; k < cnt; k++)
            {
                G_L[i][j] += A[k] * (P(z[k]) * phi_dir[k][i] * phi_dir[k][j] + Q(z[k]) * phi[k][j] * phi[k][i]);
                G[i][j] += A[k] * (phi[k][i]*phi[k][j]);
            }  
            G_L[i][j] += (aa * P(a) * phi[cnt][i] * phi[cnt][j] + bb * P(b) * phi[cnt + 1][i] * phi[cnt + 1][j]);
        }
    }
    /*
        cout << "------------------------------------------------------------------------------------------------" << endl;
        show(AA, n);
        cout << "------------------------------------------------------------------------------------------------" << endl;
        cout << "Правая часть системы:" << endl;
        show_vect(d, n);
        cout << "Система выглядит так:" << endl;
        sysout(AA, d, n);

        cout << "------------------------------------------------------------------------------------------------" << endl;
        cout << "Решение системы (искомые коэффициенты):" << endl;
        show_vect(res, n);
    */
    cout << "Матрица G:" << endl;
    show_mat(G, n);
    cout << "Матрица GL:" << endl;
    show_mat(G_L, n);
    return G_L;
}
//Скалярное произведение векторов для методов обратной итерации
double scal_vect(double* a, double* b, int n)
{
    double res = 0.0;
    for (int i = 0; i < n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}
//Макисмальный элемент вектора
double Max_el_vect(double* vec, int n)
{
    double max = vec[0];
    for (int i = 1; i < n; i++)
    {
        if (fabs(vec[i]) >= fabs(max))
        {
            max = fabs(vec[i - 1]);
        }
    }
    return max;
}
//Евклидова норма вектора
double EVnorm_vec(double* vec, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += fabs(pow(vec[i], 2));
    }
    norm = sqrt(norm);
    return norm;
}
//Апостериорная оценка собственного числа для скалярного метода 
double Apost(double** A, int n, double* Y, double lambda)
{
    double* lY = new double[n];
    for (int i = 0; i < n; i++)
    {
        A[i][i] -= lambda;
    }
    double* sup1 = Mat_x_vect(A, Y, n);

    double apost = EVnorm_vec(sup1, n) / EVnorm_vec(Y, n);
    return apost;
}
//Скалярный метод
double Scal_method(double** A, int n, double EPS)
{
    int counter = 0;
    double upper = 0.0;
    double benif = 0.0;
    double norminator = 1;
    double lambda_k = 1.0;
    double lambda_k_1 = 0.0;
    double lambda_k_2 = 0.1;
    double sup1 = 0.0;
    double sup2 = 0.0;
    double Aitken = 0.0;
    double apost = 0.0;
    double* Y = new double[2*n];
    double* Y_next = new double[2*n];
    double** H = new double* [2*n];
    for (int i = 0; i < n; i++)
    {
        H[i] = new double[n];
        Y[i] = 1;
    }
    for (int i = 0; i < n; i++)
    {
        Y[i] = 1;
    }
    std::cout << "Метод скалярных произведений." << endl;
    std::cout << "Начальный вектор Y_0: " << endl;
    show_vect(Y, n);
    cop_mat(H, A, n);
    std::cout << "|----------------|---------------|---------------|---------------|---------------|---------------|\n";
    std::cout << "|        k       |      l_k      |    l_k-l_k-1  |     apost     |     Aitken    |   Aitken-l_k  |\n";
    std::cout << "|----------------|---------------|---------------|---------------|---------------|---------------|\n";
    while (fabs(lambda_k - lambda_k_1) >= EPS)
    {
        counter++;
        cop_mat(H, A, n);
        norminator = Max_el_vect(Y, n);
        for (int i = 0; i < n; i++)
        {
            Y[i] = Y[i] / norminator;
        }
        benif = scal_vect(Y, Y, n);
        Y_next = Mat_x_vect(H, Y, n);
        upper = scal_vect(Y_next, Y, n);
        cop_vect(Y, Y_next, n);
        lambda_k_2 = lambda_k_1;
        lambda_k_1 = lambda_k;
        lambda_k = upper / benif;
        apost = Apost(H, n, Y, lambda_k);
        sup1 = (lambda_k * lambda_k_2 - pow(lambda_k_1, 2));
        sup2 = (lambda_k - 2 * lambda_k_1 + lambda_k_2);
        Aitken = sup1 / sup2;
        std::cout << "|       " << counter << "             " << lambda_k << "        " << fabs(lambda_k - lambda_k_1) << "       " << apost << "       " << Aitken << "        " << fabs(Aitken - lambda_k) << "   | \n";
        std::cout << "|------------------------------------------------------------------------------------------------|\n";
    }
    norminator = Max_el_vect(Y, n);
    for (int i = 0; i < n; i++)
    {
        Y[i] = Y[i] / norminator;
    }
    std::cout << "Соотвествующий собственный вектор: " << endl;
    show_vect(Y, n);
    std::cout << "Всего итераций: " << counter << endl;
    return lambda_k;
}
//Поиск противоположной границы спектра 
double Opposite_side(double** A, double lambda, int n, double EPS)
{
    double** B = new double* [n];
    double lambda_B = 0.0;
    for (int i = 0; i < n; i++)
    {
        B[i] = new double[n];
    }
    cop_mat(B, A, n);
    for (int i = 0; i < n; i++)
    {
        B[i][i] -= lambda;
    }
    std::cout << "Ищем противоположную границу спектра." << endl;
    std::cout << "Матрица B = A - lambda(A)*E: " << endl;
    show_mat(B, n);
    lambda_B = Scal_method (B, n, EPS);
    lambda += lambda_B;
    std::cout << "Противоположная граница спектра (Метод скалярных произведений):   " << lambda << endl;
    return lambda;
}
int main()
{
    double** a, ** H;
    
    system("chcp 1251");
    system("cls");
   // double* test1 = Legendre_coord_sys(7, 1);
   // double* test2 = Legendre_Dir_coord_sys(7, 1);
   // double* test3 = Legendre_Sec_Dir_coord_sys(7, 1);
    int n = 7;
        double** GL = Ritz_G_L(-1, 1, n);
        double* eigenvalues = Jacob(GL, n);
        double Scal = Scal_method(GL, n, 0.000001);
        double Opposite = Opposite_side(GL, Scal, n, 0.00001);
        cout << "|----------------|---------------|---------------|---------------|\n";
        cout << "|        n       | lambda_min_R  |   lambda_min  |     невязка   |\n";
        cout << "|----------------|---------------|---------------|---------------|\n";
        cout << fixed << setprecision(6) << "|        " << n << "       |   " << eigenvalues[0] << "    |    " << Opposite << "   |      " << fabs(eigenvalues[0] - Opposite) << " | \n";
        cout << "|----------------|---------------|---------------|---------------|\n";
        cout << "|----------------|---------------|---------------|---------------|\n";
        cout << "|        n       | lambda_max_R  |   lambda_max  |     невязка   |\n";
        cout << "|----------------|---------------|---------------|---------------|\n";
        cout << "|        " << n << "       |   " << eigenvalues[n - 1] << "   |    " << Scal << "  |      " << fabs(eigenvalues[n - 1] - Scal) << " | \n";
        cout << "|----------------|---------------|---------------|---------------|\n";
    
    return 0;
}