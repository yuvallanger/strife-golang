import scipy as sp
import numpy as np
import scipy.signal

def diffuse(b,c,direction):
    row = (sp.array((0,0,1,1))+c[0])%b.shape[0]
    col = (sp.array((0,1,0,1))+c[1])%b.shape[1]
    if direction:
        b[row,col] = b[row[[1,2,3,0]], col[[1,2,3,0]]]
    else:
        b[row,col] = b[row[[2,3,0,1]], col[[2,3,0,1]]]

def competiroll(N):
    # We'll use relative positions to compute exact positions of 2nd competitor cell
    NEIGHBOUR_ROW = sp.array([-1,  0,  1, -1,  0,  1, -1,  1])
    NEIGHBOUR_COL = sp.array([-1, -1, -1,  1,  1,  1,  0,  0])
    NEIGHBOUR_REL_POS = sp.array(zip(NEIGHBOUR_ROW, NEIGHBOUR_COL))
    c1 = sp.random.randint(N, size=2)
    c2 = c1 + NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)][0]
    return c1, c2

# Board size
N = 20

# radius
S_rad = 1
C_rad = 1

# diameter
dia = lambda x: 2 * x + 1
S_len = dia(S_rad)
C_len = dia(C_rad)

# the convolution matrix used to count neighbours
S_counter = sp.ones((S_len, S_len))
C_counter = sp.ones((C_len, C_len)) # convolution matrix used to count cells that produce public goods

# neighbours threshold
S_th = 3 # quorum threshold
C_th = 3 # Cooperation threshold. Above it, public goods makes a difference.

# A cell can be Signalling and/or Receptive and/or Cooperative
S = sp.rand(N, N) < 0.5
R = sp.rand(N, N) < 0.5
C = sp.rand(N, N) < 0.5

# we'll increase this by one every time two cells compete.
tick = 0

while True:
    competitor_1, competitor_2 = competiroll(N)
    # competitor_2's coordinates in a torus
    competitor_2t = competitor_2 % N 
    while ((R[competitor_1[0],competitor_1[1]] == R[competitor_2t[0], competitor_2t[1]]) and
           (S[competitor_1[0],competitor_1[1]] == S[competitor_2t[0], competitor_2t[1]]) and
           (C[competitor_1[0],competitor_1[1]] == C[competitor_2t[0], competitor_2t[1]])):
        # we'll run this until we get a new 
        competitor_1, competitor_2 = competiroll(N)
        competitor_2t = competitor_2 % N
        tick += 1
    print "competitor_1, competitor_2"
    print competitor_1, competitor_2
    # rl: row low; rh: row high
    # cl: col low; ch: col high
    rl, rh = sp.sort((competitor_1[0], competitor_2[0]))
    cl, ch = sp.sort((competitor_1[0], competitor_2[0]))
    ## produce torusified versions of boards:
    # sub array of Signallers board:
    S_sub = S[sp.arange(rl - S_rad - C_rad, rh + S_rad + C_rad)%N,:][:,sp.arange(cl - S_rad - C_rad, ch + S_rad + C_rad)%N]
    R_sub = R[sp.arange(rl - C_rad, rh + C_rad)%N, :][:, sp.arange(cl - C_rad, ch + C_rad)%N]
    C_sub = C[sp.arange(rl - C_rad, rh + C_rad)%N, :][:, sp.arange(cl - C_rad, ch + C_rad)%N]
    print "S_sub.shape, R_sub.shape, C_sub.shape"
    print S_sub.shape, R_sub.shape, C_sub.shape
    # per cell signal available
    S_conv = sp.signal.convolve2d(S_sub, S_counter, mode='valid')  
    cooping_nh = (C_sub == R_sub) == (S_conv > S_th)
    # how many cooperators around each competitor?
    print "cooping_nh"
    print cooping_nh.shape
    print cooping_nh
    C_conv = sp.signal.convolve2d(cooping_nh, C_counter, mode='valid')
    tick += 1
