import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# File name
filename = 'inputs/large_plot.txt'
COLORMAP_NAME = 'cool'

# Read and parse the data
data_lines = []
with open(filename, 'r') as f:
    lines = f.readlines()
    # Skip the first line (header)
    for i, line in enumerate(lines):
        if i == 0:
            continue
        parts = line.split()
        # Skip the first column, which is the thread number
        numeric_values = [p.replace(',', '.') for p in parts[1:]]
        numeric_floats = list(map(float, numeric_values))
        data_lines.append(numeric_floats)

data = np.array(data_lines)

# Determine shape
num_threads, num_procs = data.shape

# Create the coordinate arrays
threads = np.arange(1, num_threads + 1)
procs = np.arange(1, num_procs + 1)
X, Y = np.meshgrid(procs, threads)

xpos = X.flatten()
ypos = Y.flatten()
zpos = np.zeros_like(xpos)
dz = data.flatten()

# Reduce bar width to avoid overlap
dx = 0.6 * np.ones_like(zpos)
dy = 0.6 * np.ones_like(zpos)

# Normalize data for colormap
norm = plt.Normalize(dz.min(), dz.max())

# Use a colormap that transitions from red to dark blue
colors = plt.colormaps[COLORMAP_NAME](norm(dz))

# Create figure and 3D axes
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D bars
ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, shade=True)

# Set labels
ax.set_xlabel('Number of Processors')
ax.set_ylabel('Number of Threads')
ax.set_zlabel('Execution Time (seconds)')

ax.set_xticks(procs)
ax.set_yticks(threads)

# Set z-axis ticks from 0.* intervals
z_ticks = np.arange(0, dz.max() + 0.75, 0.75)
ax.set_zticks(z_ticks)

# Add horizontal dotted gridlines
ax.zaxis.grid(True, linestyle='--')

# Create the colorbar and associate it with the current axes
m = cm.ScalarMappable(cmap=plt.colormaps[COLORMAP_NAME], norm=norm)
m.set_array(dz)
fig.colorbar(m, ax=ax, shrink=0.5, aspect=10, label='Execution Time')

# Adjust view angle if desired
ax.view_init(elev=30, azim=65)

plt.title("Execution Time as a function of Processors and Threads (Large Dataset = 1100)")
plt.tight_layout()
plt.savefig('../plots_png/QR_plot_times_cores_hybrid_large.png')