#include "task.h"

using namespace std;

int main() {

    int N;
    int det = 0;

    srand(time(NULL));

    cout << "Enter N: ";
    cin >> N;

    double **matrix = new double *[N];

    for (int i = 0; i < N; i++)
        matrix[i] = new double [N];

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            //cout << "Enter matrix[" << i << "][" << j << "] = ";
            cin >> matrix[i][j];
        }
    }

    /*for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            matrix[i][j] = GetRandomNumber(1,10);
        }
    }*/

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }

    det = Determinant(matrix, N);
    cout << "Determinant = " << det << endl;

    if(det == 0.){
        cout << "The inverse matrix does not exist because the determinant is 0!";
    }
    else{
        Inversion(matrix, N);
        cout << "Inverse matrix:" << endl;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++)
                cout << matrix[i][j] << "  ";
            cout << endl;
        }
    }

    for (int i = 0; i < N; i++)
        delete [] matrix[i];

    delete [] matrix;

    cin.get();
    return 0;
}

