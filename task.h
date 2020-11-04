#include <iostream>
#include <stdlib.h>
#include <time.h>

int GetRandomNumber(int min, int max); // Генерация случайного числа в диапазоне [min, max]

void PrintMatrix(double **Matrix, int n); // Вывод матрицы

void GetMatrix(double **matrix, int size, int row, int col, double **newMatrix); //Матрица matrix без row-ой строки и col-того столбца, результат в newMatrix

double Determinant(double **matrix, int size); // Определитель матрицы

void Inversion(double **A, int N); // Нахождение обратной матрицы
