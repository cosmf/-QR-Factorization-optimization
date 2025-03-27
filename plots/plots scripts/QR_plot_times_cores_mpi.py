import matplotlib.pyplot as plt

# Data
cores = [1, 2, 4, 8]
serial_time_small = [8.870] * len(cores)        #800
serial_time_medium = [21.850] * len(cores)      #900
serial_time_large = [25.659] * len(cores)       #1100

mpi_times_small = [1.740, 1.230, 0.791, 0.991]      #800
mpi_times_medium = [2.440, 1.606, 1.097, 1.126]     #900
mpi_times_large = [4.652, 3.486, 2.296, 1.704]     #1100

# Calculate speedup
def calculate_speedup(serial_time, parallel_times):
    return [serial_time / p for p in parallel_times]

# Calculate parallelizable portion (P) using Amdahl's Law
def calculate_parallel_fraction(speedups, cores):
    fractions = []
    for core, speedup in zip(cores[1:], speedups[1:]):  # Skip 1 core (no parallelism)
        P = (core * (speedup - 1)) / (speedup * (core - 1))
        fractions.append(P)
    return fractions

# Compute speedups
mpi_speedups_small = calculate_speedup(serial_time_small[0], mpi_times_small)
mpi_speedups_medium = calculate_speedup(serial_time_medium[0], mpi_times_medium)
mpi_speedups_large = calculate_speedup(serial_time_large[0], mpi_times_large)

# Compute parallelizable portions
mpi_parallel_small = calculate_parallel_fraction(mpi_speedups_small, cores)
mpi_parallel_medium = calculate_parallel_fraction(mpi_speedups_medium, cores)
mpi_parallel_large = calculate_parallel_fraction(mpi_speedups_large, cores)

# Average parallelizable fraction
mpi_average_parallel = (
    sum(mpi_parallel_small + mpi_parallel_medium + mpi_parallel_large)
    / (3 * len(mpi_parallel_small))
)

# Print the parallelizable portions
print(f"MPI Parallelizable Portion (Small): {[f'{p:.2%}' for p in mpi_parallel_small]}")
print(f"MPI Parallelizable Portion (Medium): {[f'{p:.2%}' for p in mpi_parallel_medium]}")
print(f"MPI Parallelizable Portion (Large): {[f'{p:.2%}' for p in mpi_parallel_large]}")

# Print the results
print(f"Average parallelizable portion (MPI): {mpi_average_parallel:.2%}")

# Plot Acceleration
plt.figure(figsize=(12, 8))

# MPI Acceleration
plt.plot(cores, mpi_speedups_small, marker='o', linestyle='-', label='MPI (Small)', color='navy')
plt.plot(cores, mpi_speedups_medium, marker='o', linestyle='-', label='MPI (Medium)', color='firebrick')
plt.plot(cores, mpi_speedups_large, marker='o', linestyle='-', label='MPI (Large)', color='darkorange')

plt.title('Acceleration vs Number of Cores (MPI)', fontsize=14)
plt.xlabel('Cores', fontsize=12)
plt.ylabel('Acceleration (Speedup)', fontsize=12)
plt.axhline(y=1, color='black', linestyle='--', alpha=0.7, linewidth=1)  # Line for no speedup
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.savefig("../plots_png/acceleration_vs_cores_mpi.png", dpi=300)
plt.show()