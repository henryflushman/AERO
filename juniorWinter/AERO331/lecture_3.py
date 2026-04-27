import numpy as np
import matplotlib.pyplot as plt

from trusslib_v1_1 import *
from cross_truss import *


# ---------- trusslib ---------- #

# -- Create Cross Truss example -- #

truss3 = setup_truss('cross_truss_3.txt')
truss3.draw()
truss3.solve()
truss3.plot(magnification=10, structure=False)
truss3.print_displacements()
ratio19 = truss3.__getattribute__('stresses')[0]/truss3.__getattribute__('stresses')[8]
print(ratio19)


ASPECT_RATIO = 5
N_LATTICE = 3
create_cross_truss(ASPECT_RATIO*1.0, 1.0, N_LATTICE, 7e10, 5e-4, -1e3, 'lattice_truss.txt')

ltruss = setup_truss('lattice_truss.txt')
ltruss.solve()
ltruss.plot(magnification=100, structure=False)
ltruss.print_displacements(ASPECT_RATIO*N_LATTICE)

# -- Cross Truss example -- #

ctruss = setup_truss('C_truss.txt')
ctruss.draw()
ctruss.solve()
ctruss.plot(magnification=10, structure=False)

# -- Open Box Truss example -- #

openbox_truss = setup_truss('open_box_truss.txt')
openbox_truss.draw()
openbox_truss.solve()
openbox_truss.plot(magnification=10, structure=False)

# -- Closed Box Truss example -- #

closebox_truss = setup_truss('closed_box_truss.txt')
closebox_truss.draw()
closebox_truss.solve()
closebox_truss.plot(magnification=10, structure=False)