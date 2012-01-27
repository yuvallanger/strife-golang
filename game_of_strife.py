import scipy as sp
import numpy as np
import scipy.signal

def torusify(b, rl, rh, cl, ch):
    # rl, rh - row low, row high
    # cl, ch - col low, col high
    row_num, col_num = b.shape
    ri = arange(-rl, rh) % row_num
    ci = arange(-cl, ch) % col_num
    return b[ri,:][:,ci]

def diffuse(b):
    

# Board size
N = 20

# Signal
S_rad = 3 # radius
S_len = 2 * S_rad + 1 # one dimension's length
S_nh = ones((S_len, S_len)) # the convolution matrix used to count signalling neighbours

# Receptor
R_th = 3 # quorum threshhold

# Cooperation
C_rad = 3 # radius
C_len = 2 * C_rad + 1 # one dimension's length
C_nh = ones((C_len, C_len)) # convolution matrix used to count public goods producing cells
C_th = 3 # threshhold of public goods effect

# A cell can be Signalling and/or Receptive and/or Cooperative
R = sp.rand(N, N) < 0.5
S = sp.rand(N, N) < 0.5
C = sp.rand(N, N) < 0.5

# We'll use relative positions to compute exact positions of 2nd competitor cell
NEIGHBOUR_ROW = sp.array(2 * [-1, 0, 1] + [-1, 1])
NEIGHBOUR_COL = sp.array(6 * [-1] + [0, 0])
NEIGHBOUR_REL_POS = sp.array(zip(NEIGHBOUR_ROW, NEIGHBOUR_COL))

tick = 0

while True:
    competitor_1 = sp.random.randint(N, size=2)
    competitor_2 = competitor_1 + NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)]
    while (R[competitor_1] == R[competitor_2] and
           S[competitor_1] == S[competitor_2] and
           C[competitor_1] == C[competitor_2]):
        competitor_1 = sp.random.randint(N, size=2)
        competitor_2 = competitor_1 + NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)]
        tick += 1
    if competitor_1[0] > competitor_2[0]:
        rl, rh = competitor_2[0], competitor_1[0]
    else:
        rl, rh = competitor_1[0], competitor_2[0]
    if competitor_1[1] > competitor_2[1]:
        cl, ch = competitor_2[1], competitor_1[1]
    else:
        cl, ch = competitor_1[1], competitor_2[1]
    S_sub = torify(S, rl-S_rad-R_rad, rh+S_rad+R_rad, cl-S_rad-R_rad, ch+S_rad+R_rad) # sub array of Signallers board
    R_sub = torify(R, rl-R_rad, rh+S_rad, cl-R_rad, ch+R_rad)
    S_conv = sp.signal.convolve2d(S_nh, S_sum, mode='valid') # per cell signal available
    (R and S_conv < S_th)
