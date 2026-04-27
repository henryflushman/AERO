import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm

a = 1.0
b = 0.5
t = 8e-3

E = 7e10
nu = 0.3

P = 1e4

D = E*t**3/(12*(1-nu**2))

x = sp.symbols('x')
y = sp.symbols('y')

def a_del_Pi(w, del_w):
    w_x = sp.diff(w, x)
    w_y = sp.diff(w, y)
    
    w_xx = sp.diff(w_x, x)
    w_yy = sp.diff(w_y, y)
    w_xy = sp.diff(w_x, y)
    
    del_w_x = sp.diff(del_w, x)
    del_w_y = sp.diff(del_w, y)
    
    del_w_xx = sp.diff(del_w_x, x)
    del_w_yy = sp.diff(del_w_y, y)
    del_w_xy = sp.diff(del_w_x, y)
    
    a_term_1 = D * (w_xx + w_yy) * (del_w_xx + del_w_yy)
    a_term_2 = D * (1 - nu) * (
        w_xx * del_w_yy + w_yy * del_w_xx
        - 
        2 * w_xy * del_w_xy
    )
    
    return sp.integrate(
        sp.integrate(
            a_term_1 - a_term_2, (y, -b/2, b/2)
        ),
        (x, -a/2, a/2)
    )
    
def l_del_Pi(del_w):
    integrand = P * del_w
    return sp.integrate(
        sp.integrate(
            integrand, (y, -b/2, b/2)
        ),
        (x, -a/2, a/2)
    )

phi = [
    (x**2 - a**2/4)*(y**2 - b**2/4),
    (x**2 - a**2/4)*(y**2 - b**2/4)*x**2,
    (x**2 - a**2/4)*(y**2 - b**2/4)*x**4,
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**2,
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**4,
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**2*x**2,
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**4*x**4,
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**2*x**2*(x**2 + y**2),
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**2*x**2*(x**4 + y**4),
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**4*x**4*(x**2 + y**2),
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**4*x**2*(x**4 + y**4),
    (x**2 - a**2/4)*(y**2 - b**2/4)*y**4*x**2*(x**4 + y**4)*(x**2 + y**2)
]
n_basis = len(phi)

K = np.zeros((n_basis, n_basis))

for i in range(n_basis):
    for j in range(i, n_basis):
        K[i, j] = a_del_Pi(phi[j], phi[i])

for i in range(n_basis):
    for j in range(i):
        K[i, j] = K[j, i]
        
F = np.zeros(n_basis)

for i in range(n_basis):
    F[i] = l_del_Pi(phi[i])

c = np.linalg.solve(K, F)
print(f'Coefficients: ')
for i in range(n_basis):
    print(c[i])
    
    
w = c[0] * phi[0]    
for i in range(1, n_basis):
    w += c[i] * phi[i]

f_w = sp.lambdify((x, y), w, "numpy")

print(f"Midplate deflection = {1000*f_w(0.0, 0.0):.3f} mm")

xs = np.linspace(-a/2, a/2, 20)
ys = np.linspace(-b/2, b/2, 20)
xs, ys = np.meshgrid(xs, ys)

ws = f_w(xs, ys)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(xs, ys, 100*ws, cmap=cm.viridis, lw=0)
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('w (mm)')
ax.set_box_aspect((np.ptp(xs), np.ptp(ys), np.ptp(100*ws)))
plt.show()