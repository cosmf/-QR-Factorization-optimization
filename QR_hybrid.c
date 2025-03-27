#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h> // Include OpenMP header

// Function to calculate the norm of a vector using OpenMP
double vector_norm(double *v, int n) {
    double norm = 0.0;
    #pragma omp parallel for reduction(+:norm)
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}

// Function to perform QR Factorization using Householder Reflectors in parallel
void qr_factorization(int N, double *A_local, int rows_per_proc, int rank, int size) {
    double *v = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));

    // Temporary buffer to gather the k-th column
    double *column_k = (double *)malloc(N * sizeof(double));

    for (int k = 0; k < N - 1; k++) {
        // Step 1: Each process extracts its portion of the k-th column
        double local_element = 0.0;
        for (int i = 0; i < rows_per_proc; i++) {
            int global_i = rank * rows_per_proc + i;
            if (global_i == k) {
                local_element = A_local[i * N + k];
                break; // Found the element, no need to continue
            }
        }

        // Step 2: Gather the entire k-th column to rank 0
        MPI_Gather(&local_element, 1, MPI_DOUBLE, column_k, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double norm_v;
        // Step 3: Rank 0 constructs the Householder vector
        if (rank == 0) {
            // Initialize v
            #pragma omp parallel for
            for (int i = 0; i < N; i++) {
                v[i] = column_k[i];
            }

            // Compute the norm of the vector v[k:N]
            double norm = vector_norm(&v[k], N - k);
            v[k] += (v[k] >= 0) ? norm : -norm;

            // Normalize the Householder vector u
            norm_v = vector_norm(&v[k], N - k);
            
            #pragma omp parallel for
            for (int i = 0; i < N; i++) {
                u[i] = (i < k) ? 0.0 : v[i] / norm_v;
            }
        }

        // Step 4: Broadcast the Householder vector u to all processes
        MPI_Bcast(&u[k], N - k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Step 5: Each process applies the reflector to its local rows
        // Only modify rows >= k
        // We have an outer loop over i (rows) and an inner loop over j (columns) -> use collapse(2) if beneficial
        #pragma omp parallel for
        for (int i = 0; i < rows_per_proc; i++) {
            int global_i = rank * rows_per_proc + i;
            if (global_i >= k) {
                double dot = 0.0;
                for (int j = k; j < N; j++) {
                    dot += u[j] * A_local[i * N + j];
                }
                for (int j = k; j < N; j++) {
                    A_local[i * N + j] -= 2.0 * u[j] * dot;
                }
            }
        }

        // Synchronize all processes
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(v);
    free(u);
    free(column_k);
}

int main(int argc, char *argv[]) {
    int rank, size;
    int N;

    // Initialize MPI with thread support
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        printf("The MPI library does not have enough thread support\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set the number of OpenMP threads from command line if provided
    int num_threads = 8; // default number of threads
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }
    omp_set_num_threads(num_threads);

    double start_time, end_time;

    if (rank == 0) {
        // Only process 0 reads the matrix size
        N = 800;
    }

    // Broadcast the matrix size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Determine the number of rows per process
    int rows_per_proc = N / size;
    int remainder = N % size;
    if (rank < remainder) {
        rows_per_proc += 1;
    }

    // Allocate memory for the local portion of A
    double *A_local = (double *)malloc(rows_per_proc * N * sizeof(double));
    if (A_local == NULL) {
        fprintf(stderr, "Memory allocation failed for A_local on rank %d.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize matrix A on rank 0
    double *A = NULL;
    if (rank == 0) {
        A = (double *)malloc(N * N * sizeof(double));
        if (A == NULL) {
            fprintf(stderr, "Memory allocation failed for A on rank 0.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        srand((unsigned int)time(NULL));
        // Parallel initialization of A
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i * N + j] = (double)(rand() % 10 + 1);
            }
        }
    }

    // Scatter the rows of A to all processes
    // First, calculate send counts and displacements
    int *sendcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        sendcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        if (sendcounts == NULL || displs == NULL) {
            fprintf(stderr, "Memory allocation failed for sendcounts/displs on rank 0.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 0; i < size; i++) {
            sendcounts[i] = (N / size) * N;
            if (i < (N % size)) {
                sendcounts[i] += N; // Add one extra row
            }
            displs[i] = 0;
            for (int j = 0; j < i; j++) {
                displs[i] += sendcounts[j];
            }
        }
    }

    // Scatter the data
    MPI_Scatterv(
        A,
        sendcounts,
        displs,
        MPI_DOUBLE,
        A_local,
        rows_per_proc * N,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    // Start timing the QR factorization
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Perform QR factorization in parallel (MPI + OpenMP)
    qr_factorization(N, A_local, rows_per_proc, rank, size);

    // End timing
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    double time_qr = end_time - start_time;

    // Gather the results back to the root process
    MPI_Gatherv(
        A_local,
        rows_per_proc * N,
        MPI_DOUBLE,
        A,
        sendcounts,
        displs,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    // Only the root process prints the timing
    if (rank == 0) {
        printf("Time for QR factorization: %f seconds\n", time_qr);
    }

    // Free allocated memory
    free(A_local);
    if (rank == 0) {
        free(A);
        free(sendcounts);
        free(displs);
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
