#include "task.h"

using namespace std;

int GetRandomNumber(int min, int max) {
    int num = min + rand() % (max - min + 1);
    return num;
}


void PrintMatrix(double **Matrix, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++){
            cout << Matrix[i][j] << " ";
        }
        cout << endl;
    }
}


void GetMatrix(double **matrix, int size, int row, int col, double **newMatrix) {
    int offsetRow = 0; //�������� ������� ������ � �������
    int offsetCol = 0; //�������� ������� ������� � �������
    for(int i = 0; i < size-1; i++) {
        //���������� row-�� ������
        if(i == row) {
            offsetRow = 1; //��� ������ ��������� ������, ������� ���� ����������, ������ �������� ��� �������� �������
        }

        offsetCol = 0; //�������� �������� �������
        for(int j = 0; j < size-1; j++) {
            //���������� col-�� �������
            if(j == col) {
                offsetCol = 1; //��������� ������ �������, ��������� ��� ���������
            }

            newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
        }
    }
}


double Determinant(double **matrix, int size) {
    double det = 0;
    int degree = 1; // (-1)^(1+j) �� ������� ������������

    //������� ������ �� ��������
    if(size == 1) {
        return matrix[0][0];
    }
    //������� ������ �� ��������
    else if(size == 2) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else {
        //������� ��� ������ � �������
        double **newMatrix = new double*[size-1];
        for(int i = 0; i < size-1; i++) {
            newMatrix[i] = new double[size-1];
        }

        //������������ �� 0-�� ������, ���� ����� �� ��������
        for(int j = 0; j < size; j++) {
            //������� �� ������� i-� ������ � j-�� �������
            //��������� � newMatrix
            GetMatrix(matrix, size, 0, j, newMatrix);

            //����������� �����
            //�� �������: ����� �� j, (-1)^(1+j) * matrix[0][j] * minor_j (��� � ���� ����� �� �������)
            //��� minor_j - �������������� ����� �������� matrix[0][j]
            // (�������, ��� ����� ��� ������������ ������� ��� 0-�� ������ � j-�� �������)
            det = det + (degree * matrix[0][j] * Determinant(newMatrix, size-1));
            //"�����������" ������� ���������
            degree = -degree;
        }

        //������ ������ �� ������ ���� ��������(�����!)
        for(int i = 0; i < size-1; i++) {
            delete [] newMatrix[i];
        }
        delete [] newMatrix;
    }

    return det;
}


void Inversion(double **A, int N) {

double temp;

    double **E = new double *[N];

    for (int i = 0; i < N; i++)
        E[i] = new double [N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++){
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++){
        temp = A[k][k];

        for (int j = 0; j < N; j++){
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++){
            temp = A[i][k];

            for (int j = 0; j < N; j++){
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--){
        for (int i = k - 1; i >= 0; i--){
            temp = A[i][k];

            for (int j = 0; j < N; j++){
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}


