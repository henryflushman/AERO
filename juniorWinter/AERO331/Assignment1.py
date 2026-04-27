import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 2, 1001)
V = -5*x**2 + 20*x - 20
M = -(5/3)*x**3 + 10*x**2 - 20*x + (40/3)

plt.figure()
plt.plot(x, V, label="Shear Force V(x)")
plt.xlim([0, 2])
plt.xlabel("x (m)")
plt.ylim([-21, 10])
plt.ylabel("Shear Force (N)")
plt.show()

plt.figure()
plt.plot(x, M, label="Bending Moment M(x)")
plt.xlim([0, 2])
plt.xlabel("x (m)")
plt.ylim([-10, 15])
plt.ylabel("Bending Moment (N*m)")
plt.show()