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
    int offsetRow = 0; //Смещение индекса строки в матрице
    int offsetCol = 0; //Смещение индекса столбца в матрице
    for(int i = 0; i < size-1; i++) {
        //Пропустить row-ую строку
        if(i == row) {
            offsetRow = 1; //Как только встретили строку, которую надо пропустить, делаем смещение для исходной матрицы
        }

        offsetCol = 0; //Обнулить смещение столбца
        for(int j = 0; j < size-1; j++) {
            //Пропустить col-ый столбец
            if(j == col) {
                offsetCol = 1; //Встретили нужный столбец, проускаем его смещением
            }

            newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
        }
    }
}


double Determinant(double **matrix, int size) {
    double det = 0;
    int degree = 1; // (-1)^(1+j) из формулы определителя

    //Условие выхода из рекурсии
    if(size == 1) {
        return matrix[0][0];
    }
    //Условие выхода из рекурсии
    else if(size == 2) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else {
        //Матрица без строки и столбца
        double **newMatrix = new double*[size-1];
        for(int i = 0; i < size-1; i++) {
            newMatrix[i] = new double[size-1];
        }

        //Раскладываем по 0-ой строке, цикл бежит по столбцам
        for(int j = 0; j < size; j++) {
            //Удалить из матрицы i-ю строку и j-ый столбец
            //Результат в newMatrix
            GetMatrix(matrix, size, 0, j, newMatrix);

            //Рекурсивный вызов
            //По формуле: сумма по j, (-1)^(1+j) * matrix[0][j] * minor_j (это и есть сумма из формулы)
            //где minor_j - дополнительный минор элемента matrix[0][j]
            // (напомню, что минор это определитель матрицы без 0-ой строки и j-го столбца)
            det = det + (degree * matrix[0][j] * Determinant(newMatrix, size-1));
            //"Накручиваем" степень множителя
            degree = -degree;
        }

        //Чистим память на каждом шаге рекурсии(важно!)
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


