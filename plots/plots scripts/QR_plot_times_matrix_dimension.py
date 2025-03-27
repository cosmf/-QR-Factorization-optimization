import matplotlib.pyplot as plt

# Data from the file
num_points = [800, 900, 1100]
serial_times = [3.420, 5.616, 9.747]
openmp_times = [2.439, 2.599, 3.166]
mpi_times = [0.991, 1.126, 1.704]
hybrid_times = [0.700, 0.947, 1.501]

# Plotting the data with filled areas
plt.figure(figsize=(10, 6))

plt.fill_between(num_points, serial_times, alpha=0.2, label='Serial (shaded)', color='blue')
plt.fill_between(num_points, openmp_times, alpha=0.2, label='OpenMP (shaded)', color='red')
plt.fill_between(num_points, mpi_times, alpha=0.2, label='MPI (shaded)', color='green')
plt.fill_between(num_points, hybrid_times, alpha=0.2, label='Hybrid (shaded)', color='purple')

# Plotting lines on top
plt.plot(num_points, serial_times, marker='o', linestyle='-', color='blue', label='Serial')
plt.plot(num_points, openmp_times, marker='o', linestyle='-', color='red', label='OpenMP (8 threads)')
plt.plot(num_points, mpi_times, marker='o', linestyle='-', color='green', label='MPI (8 processes)')
plt.plot(num_points, hybrid_times, marker='o', linestyle='-', color='purple', label='Hybrid (MPI + OpenMP)')

# Adding titles and labels
plt.title('Execution Time vs Matrix Dimension', fontsize=14)
plt.xlabel('Matrix Dimension (NxN)', fontsize=12)
plt.ylabel('Execution Time (seconds)', fontsize=12)

# Adding grid and legend
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12, loc='upper left')

# Save the plot
plt.savefig("../plots_png/execution_time_vs_matrix_dimension_hybrid.png", dpi=300)
plt.show()