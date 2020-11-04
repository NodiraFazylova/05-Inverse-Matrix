#include <iostream>
#include <stdlib.h>
#include <time.h>

int GetRandomNumber(int min, int max); // ��������� ���������� ����� � ��������� [min, max]

void PrintMatrix(double **Matrix, int n); // ����� �������

void GetMatrix(double **matrix, int size, int row, int col, double **newMatrix); //������� matrix ��� row-�� ������ � col-���� �������, ��������� � newMatrix

double Determinant(double **matrix, int size); // ������������ �������

void Inversion(double **A, int N); // ���������� �������� �������
