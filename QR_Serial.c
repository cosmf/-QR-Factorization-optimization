#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Funcție pentru a calcula norma unui vector
double vector_norm(double *v, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}

// Funcție pentru a efectua Factorizarea QR folosind reflectori Householder
void qr_factorization(int N, double **A, double **Q, double **R) {
    double *v = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));
    double norm_v;

    // Inițializează R ca fiind A și Q ca fiind matricea identitate
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Aplicați reflectorii Householder pentru fiecare coloană a lui R
    for (int k = 0; k < N - 1; k++) {
        // Formați vectorul v
        for (int i = 0; i < N; i++) {
            v[i] = (i < k) ? 0.0 : R[i][k];
        }
        v[k] += (v[k] >= 0) ? vector_norm(v + k, N - k) : -vector_norm(v + k, N - k);
        
        // Calculați norma și formați vectorul unitate u
        norm_v = vector_norm(v + k, N - k);
        for (int i = 0; i < N; i++) {
            u[i] = (i < k) ? 0.0 : v[i] / norm_v;
        }

        // Actualizați R aplicând reflectorul Householder
        for (int i = k; i < N; i++) {
            for (int j = k; j < N; j++) {
                R[i][j] -= 2.0 * u[i] * (u[j] * R[k][j]);
            }
        }

        // Actualizați Q
        for (int i = 0; i < N; i++) {
            for (int j = k; j < N; j++) {
                Q[i][j] -= 2.0 * u[i] * (u[j] * Q[k][j]);
            }
        }
    }

    free(v);
    free(u);
}

// Funcție pentru a afișa o matrice
void print_matrix(int N, double **M, const char *name) {
    printf("Matricea %s:\n", name);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%10.4f ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main() {
    int N;
    printf("Introdu dimensiunea matricei: ");
    scanf("%d", &N);

    srand(time(NULL));

    double **A = (double **)malloc(N * sizeof(double *));
    double **Q = (double **)malloc(N * sizeof(double *));
    double **R = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
        Q[i] = (double *)malloc(N * sizeof(double));
        R[i] = (double *)malloc(N * sizeof(double));
    }

    // Generarea elementelor matricei A aleatoriu între 1 și 10
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = (double)(rand() % 10 + 1);
        }
    }

    // Afișarea matricei A
    // print_matrix(N, A, "A");

    // Efectuați factorizarea QR
    // Măsurarea timpului pentru factorizarea QR
    clock_t start_qr = clock();
    qr_factorization(N, A, Q, R);
    clock_t end_qr = clock();
    double time_qr = (double)(end_qr - start_qr) / CLOCKS_PER_SEC;
    printf("Timpul pentru factorizarea QR: %f secunde\n\n", time_qr);

    // // Măsurarea timpului pentru afișarea matricilor Q și R
    // clock_t start_print = clock();
    // clock_t end_print = clock();
    // double time_print = (double)(end_print - start_print) / CLOCKS_PER_SEC;
    // printf("Timpul pentru afișarea matricilor Q și R: %f secunde\n\n", time_print);


    // Afișați matricile Q și R
    // print_matrix(N, Q, "Q");
    // print_matrix(N, R, "R");

    // Eliberarea memoriei alocate dinamic
    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(Q[i]);
        free(R[i]);
    }
    free(A);
    free(Q);
    free(R);

    return 0;
}
