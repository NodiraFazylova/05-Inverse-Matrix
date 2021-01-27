#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
#include <future>
#include <thread>
#include <queue>
#include <condition_variable>
#include <mutex>
#include <chrono>
#include <atomic>

using namespace std;

const int Nucleus = 4; // threads
const double EPS = 1E-15;

// ******************** TIMER PART ******************** //
#ifdef __linux__
#include <unistd.h>    //usleep()
typedef std::chrono::system_clock t_clock;    //try to use high_resolution_clock on  new linux x64 computer!
#else
typedef std::chrono::high_resolution_clock t_clock;
#pragma warning(disable:4996)
#endif

std::chrono::time_point<t_clock> start_time, stop_time = start_time; char null_char = '\0';
void timer(char *title = 0, int data_size = 1) {

    stop_time = t_clock::now();

    double us = (double)chrono::duration_cast<chrono::microseconds>(stop_time - start_time).count();

    if (title)
        printf("%s time = %7lgms = %7lg MOPs\n", title, (double)us*1e-3, (double)data_size / us);

    start_time = t_clock::now();
}
// ******************** TIMER PART ******************** //

// makes columns 0
void column_zero(vector< vector<double> > &M, vector< vector<double> > &E, int pos0, int pos1, int dim, int ord);

int GetRandomNumber(int min, int max) {
    int num = min + rand() % (max - min + 1);
    return num;
}

// inverse matrix without threads
vector< vector<double> > inverse(vector< vector<double> > M){
    if (M.size() != M[0].size()){
        cout << "ERROR in inverting a non-square matrix !" << endl;
        getchar();
        return{}; // returns a null
    }

    size_t dim = M.size(); // dimension of matrix

    int i, j, ord;

    vector< vector<double> > E(dim, vector<double>(dim, 0)); // initializes output = 0
    for (i = 0; i < dim; i++){
        E[i][i] = 1.0;
    } // identity matrix

    double diag, coef;
    double *ptrM, *ptrE, *ptrM2, *ptrE2;

    for (ord = 0; ord < dim; ord++){
        int k;

        if (fabs(M[ord][ord]) < EPS){ // if that element is 0, a line that contains a non zero is added
            for (k = ord + 1; k < dim; k++){
                if (fabs(M[k][ord]) > EPS) break;
            }

            if (k >= dim)
                return{}; // error, returns null

            for (i = 0; i < dim; i++){ // added a line without 0
                M[ord][i] += M[k][i];
                E[ord][i] += E[k][i];
            }
        }

        diag = 1.0/M[ord][ord];
        ptrE = &E[ord][0];
        ptrM = &M[ord][0];

        for (i = 0; i < dim; i++){
            *ptrE++ *= diag;
            *ptrM++ *= diag;
        }

        // using the same function but without threads
        column_zero(M, E, 0, dim, dim, ord);
    } // end ord
    return E;
}

// inverse matrix with threads
vector< vector<double> > parallel_inverse(vector< vector<double> > M){

    if (M.size() != M[0].size()){
        cout << "ERROR in inverting a non-square matrix !" << endl;
        getchar();
        return{}; // returns a null
    }

    int dim = (int) M.size();
    int i, ord;

    vector< vector<double> > E(dim, vector<double>(dim, 0)); // initializes output = 0

    for (i = 0; i < dim; i++){
        E[i][i] = 1.0;
    }

    std::thread tarea[Nucleus];
    double diag;
    double *ptrM, *ptrE;

    for (ord = 0; ord < dim; ord++){
        int k;
        if (fabs(M[ord][ord]) < EPS) { // if a diagonal element is 0, it is added a column that is not 0 the diagonal element

            for (k = ord + 1; k < dim; k++){
                if (fabs(M[k][ord]) > EPS) break;
            }

            if (k >= dim)
                return{}; // error, returns null

            for (i = 0; i < dim; i++){ // looking for a line without 0, which needs to be added to make the number non-zero, in order to avoid subsequent division by 0

                M[ord][i] += M[k][i];
                E[ord][i] += E[k][i];
            }
        }

        diag = 1.0 / M[ord][ord];

        ptrE = &E[ord][0];
        ptrM = &M[ord][0];

        for (i = 0; i < dim; i++){
            *ptrE++ *= diag;
            *ptrM++ *= diag;
        }

        int pos0 = 0, N1 = dim; // initial array position
        if ((N1 < 1) || (N1 > 5000)){
            cout << "It is detected out than 1-5000 simulations points = " << N1 << "\n ABORT or press enter to continue" << endl;
            getchar();
        }
        // cout << "Initiation of " << Nucleus << " threads" << endl;
        for (int thread = 0; thread < Nucleus; thread++){

            int pos1 = (int)((thread + 1)*N1 / Nucleus); // next position
            tarea[thread] = std::thread(column_zero, std::ref(M), std::ref(E), pos0, pos1, dim, ord);
            pos0 = pos1; // next thread will work at next point
        }

        for (int thread = 0; thread < Nucleus; thread++){
            tarea[thread].join();
            // cout << "Thread number: " << thread << " end\n";
        }
    } // end ord
    return E;
}

// makes columns 0
void column_zero(vector< vector<double> > &M, vector< vector<double> > &E, int pos0, int pos1, int dim, int ord){

    double coef;
    double *ptrM, *ptrE, *ptrM2, *ptrE2;
    // we make '0' the ord column, except for the diagonal element:

    for (int i = pos0; i < pos1; i++){ // begin to end for every thread

        if (i == ord) continue;

        coef = M[i][ord]; // element to make 0

        if (fabs(coef) < EPS) continue; // if already zero, it is avoided

        ptrE = &E[i][0];
        ptrE2 = &E[ord][0];
        ptrM = &M[i][0];
        ptrM2 = &M[ord][0];

        for (int j = 0; j < dim; j++){
            *ptrE++ = *ptrE - coef * (*ptrE2++); // 1ª matrix
            *ptrM++ = *ptrM - coef * (*ptrM2++); // 2ª matrix
        }
    }
}

void TestInverseMatrix(int dim){

    cout << "|";
    cout << setw(20) << dim;
    cout << "|";
    vector< vector<double> > matrix(dim, vector<double>(dim));

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++){
            //matrix[i][j] = (-1.0 + 2.0*rand() / RAND_MAX) * 10000;
            matrix[i][j] = GetRandomNumber(0,9);
        }

    vector< vector<double> > inv_matrix, parallel_inv_matrix;

    double ini, end;

    ini = (double)clock();
    inv_matrix = inverse(matrix);
    end = (double)clock();
    cout << "|";
    cout.precision(5);
    cout << setw(25) << (end - ini) / CLOCKS_PER_SEC; // running time of the function for finding the inverse matrix without threads
    cout << "|";

    ini = end;
    parallel_inv_matrix = parallel_inverse(matrix);
    end = (double)clock();
    cout << "|";
    cout.precision(5);
    cout << setw(25) << (end - ini) / CLOCKS_PER_SEC; // running time of the function for finding the inverse matrix with thread
    cout << "||\n";

}

int main()
{
    cout << "\n   ||";
    cout << setw(15) << "Test number";
    cout << "|";

    cout << "|";
    cout << setw(20) << "Matrix dimension";
    cout << "|";

    cout << "|";
    cout << setw(25) << "Time (serial) seconds";
    cout << "|";

    cout << "|";
    cout << setw(25) << "Time (parallel) seconds";
    cout << "||\n";

    int count = 0;

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(10);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(50);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(100);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(500);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(1000);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(1500);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(2000);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(3000);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(4000);

    cout << "   ||";
    cout << setw(15) << ++count;
    cout << "|";
    TestInverseMatrix(5000);

    cout << "\n\n" << setw(75) << " ******************** The End ******************** " << endl;
    getchar();
    return 1;
}
