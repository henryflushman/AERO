import trusslib_v1_1 as truss

my_truss = truss.setup_truss('cross_truss_1.txt')

my_truss.draw()

my_truss.solve()

my_truss.plot(structure=True, magnification=1000)