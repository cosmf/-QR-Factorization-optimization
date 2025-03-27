import matplotlib.pyplot as plt

# Data
cores = [1, 2, 4, 8]
serial_time_small = [3.420] * len(cores)        #800
serial_time_medium = [5.616] * len(cores)      #900
serial_time_large = [9.747] * len(cores)       #1100

openmp_times_small = [3.496, 2.998, 3.538, 2.439]       #800
openmp_times_medium = [4.919, 3.774, 4.024, 2.599]      #900
openmp_times_large = [10.672, 5.943, 6.190, 3.166]       #1100

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
openmp_speedups_small = calculate_speedup(serial_time_small[0], openmp_times_small)
openmp_speedups_medium = calculate_speedup(serial_time_medium[0], openmp_times_medium)
openmp_speedups_large = calculate_speedup(serial_time_large[0], openmp_times_large)

# Compute parallelizable portions
openmp_parallel_small = calculate_parallel_fraction(openmp_speedups_small, cores)
openmp_parallel_medium = calculate_parallel_fraction(openmp_speedups_medium, cores)
openmp_parallel_large = calculate_parallel_fraction(openmp_speedups_large, cores)

# Average parallelizable fraction
openmp_average_parallel = (
    sum(openmp_parallel_small + openmp_parallel_medium + openmp_parallel_large)
    / (3 * len(openmp_parallel_small))
)

# Print the parallelizable portions
print(f"OpenMP Parallelizable Portion (Small): {[f'{p:.2%}' for p in openmp_parallel_small]}")
print(f"OpenMP Parallelizable Portion (Medium): {[f'{p:.2%}' for p in openmp_parallel_medium]}")
print(f"OpenMP Parallelizable Portion (Large): {[f'{p:.2%}' for p in openmp_parallel_large]}")

# Print the results
print(f"Average parallelizable portion (OpenMP): {openmp_average_parallel:.2%}")

# Plot Acceleration
plt.figure(figsize=(12, 8))

# OpenMP Acceleration
plt.plot(cores, openmp_speedups_small, marker='o', linestyle='-', label='OpenMP (Small)', color='darkblue')
plt.plot(cores, openmp_speedups_medium, marker='o', linestyle='-', label='OpenMP (Medium)', color='red')
plt.plot(cores, openmp_speedups_large, marker='o', linestyle='-', label='OpenMP (Large)', color='gold')

plt.title('Acceleration vs Number of Cores (OpenMP)', fontsize=14)
plt.xlabel('Cores', fontsize=12)
plt.ylabel('Acceleration (Speedup)', fontsize=12)
plt.axhline(y=1, color='black', linestyle='--', alpha=0.7, linewidth=1)  # Line for no speedup
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.savefig("../plots_png/acceleration_vs_cores_openmp.png", dpi=300)
plt.show()