#!/usr/bin/python

import scipy as sp
import numpy as np
import scipy.signal

## functions

def diffuse(b,c,direction):
    row = (sp.array((0,0,1,1))+c[0])%b.shape[0]
    col = (sp.array((0,1,0,1))+c[1])%b.shape[1]
    if direction:
        b[row,col] = b[row[[1,2,3,0]], col[[1,2,3,0]]]
    else:
        b[row,col] = b[row[[2,3,0,1]], col[[2,3,0,1]]]

def competiroll(N):
    """draw two competitor positions"""
    # We'll use relative positions to compute exact positions of 2nd competitor cell
    NEIGHBOUR_ROW = sp.array([-1,  0,  1, -1,  0,  1, -1,  1])
    NEIGHBOUR_COL = sp.array([-1, -1, -1,  1,  1,  1,  0,  0])
    NEIGHBOUR_REL_POS = sp.array(zip(NEIGHBOUR_ROW, NEIGHBOUR_COL))
    c1 = sp.random.randint(N, size=2)
    deltas = NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)[0]]
    c2 = c1 + deltas
    return c1, c2, deltas

## settings

# Board size
N = 20

S_cost = 3
R_cost = 8
C_cost = 30

# cooperation benefit, in ratio
benefit = 0.3

# radius
S_rad = 1
C_rad = 1

# diameter of the convolution matrix
diameter = lambda x: 2 * x + 1
S_len = diameter(S_rad)
C_len = diameter(C_rad)

# the convolution matrix used to count neighbours
S_counter = sp.ones((S_len, S_len))
C_counter = sp.ones((C_len, C_len)) # convolution matrix used to count cells that produce public goods

# neighbours effects' thresholds
S_th = 3 # quorum threshold
C_th = 3 # Cooperation threshold. Above it, public goods makes a difference.

# A cell can be Signalling and/or Receptive and/or Cooperative
S = sp.rand(N, N) < 0.5
R = sp.rand(N, N) < 0.5
C = sp.rand(N, N) < 0.5

# we'll increase this by one every time two cells compete.
tick = 0

## main stuff

while_count = 0

while True:
    print "while_count", while_count
    competitor_1, competitor_2 = competiroll(N)
    # competitor_2's coordinates in a torus:
    competitor_2t = competitor_2 % N
    # we'll run this until we get a pair of competitors that are actually different:
    while ((R[competitor_1[0],competitor_1[1]] == R[competitor_2t[0], competitor_2t[1]]) and
           (S[competitor_1[0],competitor_1[1]] == S[competitor_2t[0], competitor_2t[1]]) and
           (C[competitor_1[0],competitor_1[1]] == C[competitor_2t[0], competitor_2t[1]])):
        competitor_1, competitor_2 = competiroll(N)
        competitor_2t = competitor_2 % N
        # time passes:
        tick += 1
    print "competitor_1, competitor_2"
    print competitor_1, competitor_2
    # we figure out what are the low and high coordinates, so we may create a torusified copy of the competing cells' neighborhood.
    # "rl" stands for "row low"; "rh" for "row high"
    # "cl": "col low"; "ch": "col high"
    rl, rh = sp.sort((competitor_1[0], competitor_2[0]))
    cl, ch = sp.sort((competitor_1[1], competitor_2[1]))
    # here we produce torusified versions of the competing cells' neighborhood.
    # for signallers, we take both S_rad and C_rad around our competitors because,
    # signallers affect receptive && cooperating cells which affect our competitors
    S_sub = S[sp.arange(rl - S_rad - C_rad, rh + S_rad + C_rad + 1)%N,:][:,sp.arange(cl - S_rad - C_rad, ch + S_rad + C_rad + 1)%N]
    R_sub = R[sp.arange(rl - C_rad, rh + C_rad + 1)%N, :][:, sp.arange(cl - C_rad, ch + C_rad + 1)%N]
    C_sub = C[sp.arange(rl - C_rad, rh + C_rad + 1)%N, :][:, sp.arange(cl - C_rad, ch + C_rad + 1)%N]
    print "rl, rh, cl, ch",
    print rl, rh, cl, ch
    print "rl - S_rad - C_rad, rh + S_rad + C_rad, cl - S_rad - C_rad, ch + S_rad + C_rad"
    print rl - S_rad - C_rad, rh + S_rad + C_rad, cl - S_rad - C_rad, ch + S_rad + C_rad
    print "S_sub.shape, R_sub.shape, C_sub.shape"
    print S_sub.shape, R_sub.shape, C_sub.shape
    # we count how many signallers are within each cell's neighbourhood
    S_conv = sp.signal.convolve2d(S_sub, S_counter, mode='valid')
    # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
    cooping_cells = (C_sub == R_sub) == (S_conv > S_th)
    # how many cooperators around each competitor?
    print "cooping_cells"
    print cooping_cells.shape
    print cooping_cells
    C_conv = sp.signal.convolve2d(cooping_cells, C_counter, mode='valid')
    # Public goods effect.
    # G for Goods
    G = (C_conv > C_th)
    
    print "G.shape", G.shape
    # all cells for which the effect of goods is above threshold is True in G.
    # M for Metabolism
    if competitor_1[0] < competitor_2[0] and competitor_1[1] < competitor_2[1]:
        S_pair[competitor_1[0]:competitor_2[0]+1, competitor
    elif competitor_1[0] == competitor_2[0] and competitor_1[1] < competitor_2[1]:
    elif competitor_1[0] > competitor_2[0] and competitor_1[1] < competitor_2[1]:
    elif competitor_1[0] < competitor_2[0] and competitor_1[1] == competitor_2[1]:
    elif competitor_1[0] == competitor_2[0] and competitor_1[1] == competitor_2[1]:
    elif competitor_1[0] > competitor_2[0] and competitor_1[1] == competitor_2[1]:
    elif competitor_1[0] < competitor_2[0] and competitor_1[1] > competitor_2[1]:
    elif competitor_1[0] == competitor_2[0] and competitor_1[1] > competitor_2[1]:
    elif competitor_1[0] > competitor_2[0] and competitor_1[1] > competitor_2[1]:
    M = G * (1 - benefit) * (S_cost * S[tuple(set((competitor_1[0], + competitor_2t[0]))), :]N * competitor_2t[0] + competitor_2t[1]]] +
                             R_cost * R[[N * competitor_1[0] + competitor_1[1], N * competitor_2t[0] + competitor_2t[1]]] +
                             C_cost * C[[N * competitor_1[0] + competitor_1[1], N * competitor_2t[0] + competitor_2t[1]]])
    # all false in G don't benefit from public goods (G^True flips values)
    M += (G^True) *  (S_cost * S[[N * competitor_1[0] + competitor_1[1], N * competitor_2t[0] + competitor_2t[1]]] +
                      R_cost * R[[N * competitor_1[0] + competitor_1[1], N * competitor_2t[0] + competitor_2t[1]]] +
                      C_cost * C[[N * competitor_1[0] + competitor_1[1], N * (competitor_2t[0]) + (competitor_2t[1])]])
    tick += 1
    while_count += 1
