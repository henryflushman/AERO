import matplotlib.pyplot as plt
import numpy as np

# Debris data
pieces = ["A1", "B", "C3", "D1", "D3", "E1", "Small 1", "Small 2", "Small 3"]
distances = np.array([0.5, 1.2, 0.2, 1.7, 2.2, 2.5, 0.5, 3.2, 1.3])
angles_deg = np.array([175, 110, 85, 10, 15, 365, 210, 30, 290])

# Convert degrees to radians for polar plot
angles_rad = np.deg2rad(angles_deg)

# Create polar plot
plt.figure()
ax = plt.subplot(111, projection='polar')

# Set 0° at top (12 o'clock) and clockwise positive direction
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# Plot debris points
ax.scatter(angles_rad, distances)

# Label each fragment
for i in range(len(pieces)):
    ax.text(angles_rad[i], distances[i], pieces[i])

# Title
ax.set_title("Top-Down View of Debris Fragmentation")

plt.show()