def coords(i,j,m,n):
    coords = []
    ioffs=[0, 1]
    joffs=[1, 0]
    for ioff, joff in zip(ioffs, joffs):
        coords.append(((i+ioff)%m,(j+joff)%n))
    return coords

def step(A):
    m = A.nrows()
    n = A.ncols()
    B = copy(A)
    #netFlow = 0.1 # not a very good model.
    portion = 0.1
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            cs = coords(i,j,m,n)
            #perCellFlow = netFlow/len(cs) # lets forget that ever happened.
            for k,l in cs:
                B[i,j] -= portion * (A[i,j] - A[k,l]) # flow is proportional to difference in concentrations. Not my idea.
                B[k,l] += portion * (A[i,j] - A[k,l])
    return B

##

A = zero_matrix(RR, 10,10)
A[int(A.nrows()/3), int(A.ncols()/2)] = 1
l = []
print sum(sum(A))
a = animate([matrix_plot(A)])
for i in range(100):
    A = step(A)
    a = a * animate([matrix_plot(A)])

a.show()
