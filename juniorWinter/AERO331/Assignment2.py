import numpy as np
import matplotlib.pyplot as plt

# ----- Initialize Square ----- #

def square_boundary(n=200):
    t = np.linspace(0, 1, n)
    
    # four edges of the square
    x1, y1 = t, 0*t
    x2, y2 = 1+0*t, t
    x3, y3 = 1-t, 1+0*t
    x4, y4 = 0*t, 1-t
    
    x = np.concatenate([x1, x2, x3, x4])
    y = np.concatenate([y1, y2, y3, y4])
    return x, y

# ----- Deformations ----- #

def deform(x, y, u_func, v_func):
    return x + u_func(x, y), y + v_func(x, y)


# Situation 1
u1 = lambda x, y: 2*x - 2
v1 = lambda x, y: 3 - .25*y

# Situation 2
u2 = lambda x, y: 2*x
v2 = lambda x, y: -.25*y

# Situation 3
u3 = lambda x, y: .5*y
v3 = lambda x, y: 0*x


# ------ Simulate Deformations ------ #
x, y = square_boundary()

x1p, y1p = deform(x, y, u1, v1)
x2p, y2p = deform(x, y, u2, v2)
x3p, y3p = deform(x, y, u3, v3)

# ------ Plotting ------ #

fig, axs = plt.subplots(1, 3, figsize=(15,4))

situations = [
    ("Situation 1: u=2x-2, v=3-0.25y", x1p, y1p),
    ("Situation 2: u=2x, v=-0.25y", x2p, y2p),
    ("Situation 3: u=0.5y, v=0", x3p, y3p)
]

for ax, (title, xd, yd) in zip (axs, situations):
    # Unchanged
    ax.plot(x, y, linewidth=2, label="Original")
    # Deformed
    ax.plot(xd, yd, linewidth=2, label="Deformed")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()