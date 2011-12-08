def coords(i,j,m,n):
    """
    Will return pairs of a cell, not surrounding cells.
    """
    coords = []
    ioffs=[0, 1] # The cat's idea.
    joffs=[1, 0]
    for ioff, joff in zip(ioffs, joffs):
        coords.append(((i+ioff)%m,(j+joff)%n))
    return coords

def step(A):
    m = A.nrows()
    n = A.ncols()
    B = copy(A)
    #netFlow = 0.1 # taking % out of each cell and dividing by adjacent cells is not a good model of diffusion. A much better model is taking % out of the difference between a pair of cells, as follows.
    portion = 0.1
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            cs = coords(i,j,m,n)
            #perCellFlow = netFlow/len(cs) # lets forget that ever happened.
            for k,l in cs:
                flow = portion * (A[i,j] - A[k,l]) # flow is proportional to difference in concentrations. DMK's idea.
                B[i,j] -= flow
                B[k,l] += flow
    return B

#####
# settings
#####

steps = 100
matrixSize = (10, 20)

####
# main()-ish code
####

def diffusion_demo():
    A = zero_matrix(RR, matrixSize[0], matrixSize[1])
    A[int(A.nrows()/3), int(A.ncols()/2)] = 1
    print sum(sum(A))
    a = animate([matrix_plot(A)])
    for i in range(steps):
     A = step(A)
     a = a * animate([matrix_plot(A)])
    a.show()

diffusion_demo()
