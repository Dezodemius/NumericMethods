#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <io.h>

using namespace std;

/**
 * Тип для функций начального распределения.
 */
typedef double (*initialFuncType)(double);

/**
 * Тип для функций, которые описывают разностную схему решения ДУ.
 */
typedef void (*schemeFuncType)(int, int, double, double *const *, int);

/**
 * Функция sign. Показывает знак аргумента.
 * @param x Аргумент функции.
 * @return 0, если x < 0 и 1, если x > 0.
 */
int sgn(double x) {
    return (0.0 < x) - (x < 0.0);
}

/**
 * Явная схемы решения дифф. уравнения.
 * @param M
 * @param N
 * @param kappa
 * @param U
 */
void implicitScheme(int M, int N, double kappa, double *const *U, int nameNumber) {
    char str[19];
    sprintf(str, "data/%d_data.csv", nameNumber);
    ofstream dataFile (str);
    if (dataFile.is_open()) {
        for (int n = 0; n < M - 1; n++) {
            for (int j = 1; j < N; j++) {
                U[n + 1][j] = U[n][j] - kappa * (U[n][j] - U[n][j - 1]);
                dataFile << U[n][j] << endl;
            }
            U[n + 1][0] = U[n + 1][N - 1];
            dataFile << endl;
        }
        for (int j = 0; j < N; j++) {
            dataFile << U[M - 1][j] << endl;
        }
        dataFile.close();
    }
    else{
        cout << "Cannot open the file" << endl;
    }
}

/**
 * Схема с ограничителем Ван Лира.
 * @param M
 * @param N
 * @param kappa
 * @param U
 */
void schemeWithVanLeerLimiter(int M, int N, double kappa, double *const *U, int nameNumber) {
    char str[32];
    sprintf(str, "data/%d_data_van_leer.csv", nameNumber);
    ofstream dataFile (str);
    if (dataFile.is_open()) {
        for (int n = 0; n < M - 1; n++) {
            for (int j = 2; j < N; j++) {
                double dU1 = (U[n][j + 1] - U[n][j]) * (U[n][j] - U[n][j - 1]) <= 0.0
                        ? 0.0
                        : min(min(
                                2.0 * abs(U[n][j] - U[n][j - 1]),
                                abs(U[n][j + 1] - U[n][j - 1]) / 2.0),
                              2.0 * abs(U[n][j + 1] - U[n][j])) * sgn(U[n][j + 1] - U[n][j]);

                double dU2 = (U[n][j] - U[n][j - 1]) * (U[n][j - 1] - U[n][j - 2]) <= 0.0
                             ? 0.0
                             : min(min(
                                2.0 * abs(U[n][j - 1] - U[n][j - 2]),
                                abs(U[n][j] - U[n][j - 2]) / 2.0),
                                   2.0 * abs(U[n][j] - U[n][j - 1])) * sgn(U[n][j] - U[n][j - 1]);

                double U1 = U[n][j] + (1.0 - kappa) / 2.0 * (dU1);
                double U2 = U[n][j - 1] + (1.0 - kappa) / 2.0 * (dU2);

                U[n + 1][j] = U[n][j] - kappa * (U1 - U2);
                dataFile << U[n][j] << endl;
            }
            U[n + 1][0] = U[n + 1][N - 1];
            U[n + 1][1] = U[n + 1][N - 2];
            dataFile << endl;
        }
        for (int j = 0; j < N; j++) {
            dataFile << U[M - 1][j] << endl;
        }
        dataFile.close();
    }
    else{
        cout << "Cannot open the file" << endl;
    }
}

/**
 * Ступенчатая функция.
 * @param x Аргумент функции.
 * @param a Левая граница интервала.
 * @param b Правая граница интервала.
 * @return Значение функции в точке x.
 */
double step_function(double x, double a, double b){
    if (x >= a && x <= b)
        return 1.0;
    return 0.0;
}

/**
 * Ступенчатая функция.
 * @param x Аргумент функции.
 * @return Значение функции в точке x.
 */
double step_function(double x){
    return step_function(x, 0.2, 0.4);
}

/**
 * Треугольная функция.
 * @param x Аргумент функции.
 * @return Значение функции в точке x.
 */
double triangle_function(double x){
    if (abs(x) < 1.0)
        return 1 - abs(x);
    else
        return 0.0;
}

/**
 * Функция синуса для интервала.
 * @param x Аргумент функции.
 * @param a Левая граница интервала.
 * @param b Правая граница интервала.
 * @return Значение функции в точке x.
 */
double sin_function(double x, double a, double b){
    if (x >= a && x <= b)
        return sin(x);
    return 0.0;
}

/**
 * Функция синуса для интервала.
 * @param x Аргумент функции.
 * @return Значение функции в точке x.
 */
double sin_function(double x){
    return sin_function(x, 0.2, 0.4);
}

/**
 * Построение точного решения.
 * @param a Коэффициент уравнения.
 * @param T0 Начальное время.
 * @param T1 Конечное время.
 * @param X0 Начало координаты.
 * @param X1 Конец координаты.
 * @param h Шаг по оси X.
 * @param tau Шаг по оси T.
 * @param initialFunc Функция начального распределения.
 */
void exactSolution(double a, double T0, double T1, double X0, double X1, double h, double tau,
                   initialFuncType initialFunc, int nameNumber) {
    int M = floor((double)((T1 - T0) / tau));
    int N = floor((double)((X1 - X0) / h));
    cout << "M: " << M << endl;
    cout << "N: " << N << endl;

    auto **U1 = new double*[M];
    for (int i = 0; i < M; i++)
        U1[i] = new double[N];

    char str[25];
    sprintf(str, "data/%d_data_exact.csv", nameNumber);
    ofstream dataFile (str);
    if (dataFile.is_open()) {
        double t = T0;
        for (int n = 0; n < M; n++) {
            double x = X0;
            for (int j = 0; j < N; j++) {
                U1[n][j] = initialFunc(x - a * t);
                dataFile << U1[n][j] << endl;
                x += h;
            }
            dataFile << endl;
            U1[n][0] = U1[n][N - 1];
            t += tau;
        }
        dataFile.close();
    }
    else{
        cout << "Cannot open the file" << endl;
    }
}

/**
 * Функция для решения ДУ.
 * @param a Коэффициент в уравнении.
 * @param T0 Начальное время.
 * @param T1 Конечное время.
 * @param X0 Начало координаты.
 * @param X1 Конец координаты.
 * @param h Шаг по оси X;
 * @param tau Шаг по оси T;
 * @param schemeFunc Функция для расчёта по разностной схеме.
 * @param initialFunc Функция начального распределения.
 * @param nameNumber Номер. Нужен для записи в файл.
 */
void DiffEqSol(double a, double T0, double T1, double X0, double X1, double h, double tau, schemeFuncType schemeFunc,
               initialFuncType initialFunc, int nameNumber){
    int M = floor((double)((T1 - T0) / tau));
    int N = floor((double)((X1 - X0) / h));
    cout << "M: " << M << endl;
    cout << "N: " << N << endl;

    // Число Куранта.
    double kappa = a * tau / h;
    cout << "kappa: " << kappa << endl;
    if (kappa >= 1.0 || kappa <= 0.0){
        cout << "Значение числа Куранта не попадает в отрезок [0; 1]:" << kappa << endl;
        return;
    }

    auto **U = new double*[M];
    for (int i = 0; i < M; i++)
        U[i] = new double[N];

    double x = 0.0;
    for (int j = 0; j < N; j++) {
        U[0][j] = initialFunc(x);
        x += h;
    }

    schemeFunc(M, N, kappa, U, nameNumber);
}

/**
 * Получить величину шага сетки по крайним значениям и количеству точек.
 * @param start Начало шага.
 * @param end Конец шага.
 * @param n Количество шагов.
 * @param step Шаг.
 */
void GetStepByCount(double start, double end, double n, double &step){
    step = (end - start) / n;
}

/**
 * Входная точка в программу.
 * @return Код завершения.
 */
int main() {
    double t0 = 0.0;
    double t1 = 1.0;
    double x0 = 0.0;
    double x1 = 1.0;

    double h = 0.01;
    GetStepByCount(x0, x1, 100, h);

    double tau = 0.01;
    double a = 0.5;

    schemeFuncType schemeFuncs[2] = {schemeWithVanLeerLimiter, implicitScheme};
    initialFuncType initialFuncs[3] = {sin_function, step_function, triangle_function};

    for (auto & schemeFunc : schemeFuncs){
        for (int i = 0; i < 3; i++) {
            DiffEqSol(a, t0, t1, x0, x1, h, tau, schemeFunc, initialFuncs[i], i + 1);
            exactSolution(a, t0, t1, x0, x1, h, tau, initialFuncs[i], i + 1);
        }
    }


    return 0;
}