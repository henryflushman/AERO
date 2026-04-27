# Refined truss demo.
# Created for AERO 331 (Winter 2026)
# 
# Amuthan Ramabathiran
# Aerospace Engineering
# Cal Poly SLO
#
# January 11, 2026


def create_cross_truss(Lx, Ly, N, E, A, F, filename=None):
    # Assumption: Lx > Ly
    assert Lx > Ly, "Lx needs to be larger than Ly."

    Nx = int(Lx/Ly)*N
    Ny = 1*N

    dx = Lx/(Nx - 1)
    dy = Ly/(Ny - 1)

    # n_joints = 2*Nx + 2*Ny - 4
    # n_members = 4*Nx + 4*Ny - 10
    
    joints = []
    members = []

    # Joints
    ## Bottom edge
    for i in range(Nx):
        x = i*dx
        y = 0.0
        joints.append([x, y])

    ## Top edge
    for i in range(Nx):
        x = i*dx
        y = Ly
        joints.append([x, y])

    ## Left edge
    for i in range(Ny - 2):
        x = 0.0
        y = (i + 1)*dy 
        joints.append([x, y])

    ## Right edge
    for i in range(Ny - 2):
        x = Lx
        y = (i + 1)*dy 
        joints.append([x, y])

    # assert len(joints) == n_joints, "Missed atleast one joint."
    n_joints = len(joints)

    # Members
    ## Bottom edge
    for i in range(Nx - 1):
        members.append([i + 1, i + 2])

    ## Top edge
    for i in range(Nx - 1):
        members.append([Nx + i + 1, Nx + i + 2])

    ## Left edge
    members.append([1, 2*Nx + 1])
    for i in range(Ny - 3):
        members.append([2*Nx + i + 1, 2*Nx + i + 2])
    members.append([2*Nx + Ny - 2, Nx + 1])

    ## Right edge
    members.append([Nx, 2*Nx + Ny - 1])
    for i in range(Ny - 3):
        members.append([2*Nx + Ny + i - 1, 2*Nx + Ny + i])
    members.append([2*Nx + 2*Ny - 4, 2*Nx])

    ## Forward diagonals
    for i in range(Ny - 2):
        members.append([2*Nx + i + 1, Nx + Ny - i - 1])

    for i in range(Nx - Ny + 1):
        members.append([i + 1, Nx + Ny + i]) 

    for i in range(Ny - 2):
        members.append([Nx - Ny + 2 + i, 2*Nx + 2*Ny - 4 - i])

    ## Backward diagonals
    for i in range(Ny - 2):
        members.append([2 + i, 2*Nx + i + 1])

    for i in range(Nx - Ny + 1):
        members.append([Ny + i, Nx + 1 + i])

    for i in range(Ny - 2):
        members.append([2*Nx + Ny - 1 + i, 2*Nx - Ny + i + 2])

    # assert len(members) == n_members, "Missed atleast one member."
    n_members = len(members)

    if filename is not None:
        with open(filename, 'w') as f:
            # Block 1
            f.write('DIM\t2\n')
            f.write('---\n')

            # Block 2
            f.write('EXTRUDE NO\n')
            f.write('PATTERN NONE\n')
            f.write('LAYERS 0\n')
            f.write('DEPTH 0.0\n')
            f.write('---\n')

            # Block 3
            f.write(f'JOINTS\t{n_joints}\n')
            for i in range(n_joints):
                x, y = joints[i]
                f.write(f'{i+1}\t{x}\t{y}\n')
            f.write('---\n')

            # Block 4
            f.write(f'MEMBERS\t{n_members}\n')
            for i in range(n_members):
                n1, n2 = members[i]
                f.write(f'{i+1}\t{n1}\t{n2}\t{E}\t{A}\n')
            f.write('---\n')

            # Block 5
            f.write(f'CONSTRAINTS\t{2*Ny}\n')
            f.write('1\t1\tX\t0.0\n')
            f.write('2\t1\tY\t0.0\n')
            for i in range(Ny - 2):
                f.write(f'{2*i + 3}\t{2*Nx + 1 + i}\tX\t0.0\n')
                f.write(f'{2*i + 4}\t{2*Nx + 1 + i}\tY\t0.0\n')
            f.write(f'{2*Ny - 1}\t{Nx + 1}\tX\t0.0\n')
            f.write(f'{2*Ny}\t{Nx + 1}\tY\t0.0\n')
            f.write('---\n')

            # Block 6
            f.write(f'LOADS 1\n')
            f.write(f'1\t{Nx}\tY\t{F}\n')

    else:
        pass # Complete code to return joints, members, ...




    



