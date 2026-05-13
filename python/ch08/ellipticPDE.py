import numpy as np

def ellipticPDE(nx, ny, dx, dy, bc, f=0):
    """
    Solve 2-dimensional elliptic PDE (Python version of ellipticPDE.m)
    
    Args:
        nx, ny: number of divisions in x and y direction
        dx, dy: x and y increments
        bc: matrix of boundary conditions (4x2 or 4x3)
        f: constant for Poisson equation (default is 0 for Laplace)
    """
    
    # Initialization and validation
    bc = np.array(bc)
    a, b = bc.shape
    if a != 4:
        raise ValueError('Invalid number of boundary conditions.')
    if b < 2 or b > 3:
        raise ValueError('Invalid boundary condition.')
        
    # If bc is 4x2 and only contains Dirichlet/Neumann (types 1 or 2), pad with zeros
    if b == 2:
        bc = np.hstack((bc, np.zeros((4, 1))))
    
    nx = int(nx)
    ny = int(ny)
    x = np.arange(nx + 1) * dx
    y = np.arange(ny + 1) * dy
    dx2 = 1 / (dx**2)
    dy2 = 1 / (dy**2)
    
    # Coefficient matrix and constant vector
    n = (nx + 1) * (ny + 1)
    A = np.zeros((n, n))
    c = np.zeros(n)
    
    # Interior nodes
    # MATLAB indices 2:nx correspond to Python 1:nx
    for j in range(1, ny):
        for i in range(1, nx):
            ind = j * (nx + 1) + i
            A[ind, ind] = -2 * (dx2 + dy2)
            A[ind, ind + 1] = dx2
            A[ind, ind - 1] = dx2
            A[ind, ind + (nx + 1)] = dy2
            A[ind, ind - (nx + 1)] = dy2
            c[ind] = f

    # Boundary conditions - Helper for indexing
    # Matlab 1-based indexing to Python 0-based conversion
    
    # 1. Lower x boundary (i=0)
    for j in range(1, ny):
        ind = j * (nx + 1) + 0
        if bc[0, 0] == 1: # Dirichlet
            A[ind, ind] = 1.0
            c[ind] = bc[0, 1]
        elif bc[0, 0] in [2, 3]: # Neumann or Robbins
            A[ind, ind] = -(3 / (2 * dx) + bc[0, 2])
            A[ind, ind + 1] = 2 / dx
            A[ind, ind + 2] = -1 / (2 * dx)
            c[ind] = bc[0, 1]

    # 2. Upper x boundary (i=nx)
    for j in range(1, ny):
        ind = j * (nx + 1) + nx
        if bc[1, 0] == 1:
            A[ind, ind] = 1.0
            c[ind] = bc[1, 1]
        elif bc[1, 0] in [2, 3]:
            A[ind, ind] = (3 / (2 * dx) - bc[1, 2])
            A[ind, ind - 1] = -2 / dx
            A[ind, ind - 2] = 1 / (2 * dx)
            c[ind] = bc[1, 1]

    # 3. Lower y boundary (j=0)
    for i in range(1, nx):
        ind = 0 * (nx + 1) + i
        if bc[2, 0] == 1:
            A[ind, ind] = 1.0
            c[ind] = bc[2, 1]
        elif bc[2, 0] in [2, 3]:
            A[ind, ind] = -(3 / (2 * dy) + bc[2, 2])
            A[ind, ind + (nx + 1)] = 2 / dy
            A[ind, ind + 2 * (nx + 1)] = -1 / (2 * dy)
            c[ind] = bc[2, 1]

    # 4. Upper y boundary (j=ny)
    for i in range(1, nx):
        ind = ny * (nx + 1) + i
        if bc[3, 0] == 1:
            A[ind, ind] = 1.0
            c[ind] = bc[3, 1]
        elif bc[3, 0] in [2, 3]:
            A[ind, ind] = (3 / (2 * dy) - bc[3, 2])
            A[ind, ind - (nx + 1)] = -2 / dy
            A[ind, ind - 2 * (nx + 1)] = 1 / (2 * dy)
            c[ind] = bc[3, 1]

    # Corner nodes
    # A(1,1) -> A[0,0]
    A[0, 0] = 1; A[0, 1] = -0.5; A[0, nx + 1] = -0.5; c[0] = 0
    # A(nx+1, nx+1) -> A[nx, nx]
    A[nx, nx] = 1; A[nx, nx - 1] = -0.5; A[nx, 2 * (nx + 1) - 1] = -0.5; c[nx] = 0
    # Lower-left corner in indexing (ny*(nx+1)+1) -> A[ny*(nx+1), ny*(nx+1)]
    idx_ll = ny * (nx + 1)
    A[idx_ll, idx_ll] = 1; A[idx_ll, idx_ll + 1] = -0.5; A[idx_ll, idx_ll - (nx + 1)] = -0.5; c[idx_ll] = 0
    # Last node A(n,n) -> A[n-1, n-1]
    A[n-1, n-1] = 1; A[n-1, n-2] = -0.5; A[n-1, n - 1 - (nx + 1)] = -0.5; c[n-1] = 0

    # Solve the linear system
    u = np.linalg.solve(A, c)
    
    # Rearrange results in matrix form U(ny+1, nx+1)
    U = u.reshape((ny + 1, nx + 1))
    
    return x, y, U