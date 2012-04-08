## A model by Dr. Avigdor Eldar based on 

import scipy as sp
import scipy.signal
import pygame

# TODO: How do I figure out what's the wanted frequency of mutation and diffusion?
# TODO: Where exactly do I "global foo"? In which scopes and namespaces does "foo" lives?
# BUG: check for sameness shouldn't be on its own loop but just at the top of the algorithm's loop so each check will be followed.immediately by one diffuse call and one mutate call.
# 

## functions

def competiroll(N):
    """draw two competitor positions"""
    # We'll use relative positions to compute exact positions of 2nd competitor cell
    NEIGHBOUR_ROW = sp.array([-1,  0,  1, -1,  0,  1, -1,  1])
    NEIGHBOUR_COL = sp.array([-1, -1, -1,  1,  1,  1,  0,  0])
    NEIGHBOUR_REL_POS = sp.array(zip(NEIGHBOUR_ROW, NEIGHBOUR_COL))
    c1 = sp.random.randint(N, size=2)
    c2 = c1 + NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)[0]]
    return c1, c2


## settings

# Board size
N = 50
cell_num = N**2

## time keeping
# number of generations the simulation will run
# each generation is defined as the average number of steps for which each cell
# on the board was in a competition once since last generation (?)

# we'll increase this by one every time two cells compete.
step = 0

generations = 50
steps_final = generations * cell_num

# number of genotypes possible
genotype_num = 8

# Cost of gene expression

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

## data sampling
# we will take a frequency sample some number of times per generation
steps_per_gen = N**2
samples_per_gen = 100
samples_num = samples_per_gen * generations
steps_per_sample = sp.uint32(sp.floor(1.0 * steps_per_gen / samples_per_gen))
sample_count = 0
# We want to know the frequency of each genotype per generation
samples_frequency = sp.empty((samples_num, genotype_num))
samples_nhood = sp.empty((samples_num, genotype_num, genotype_num))

## main stuff

# pygame initialization

#pygame.init()
#screen = pygame.display.set_mode((N*4, N*4))
#pygame.display.set_caption("lets see")

def competition():
    global C, R, S, step
    ## compete
    competitor_1, competitor_2 = competiroll(N)
    # competitor_2's coordinates in a torus:
    competitor_2t = competitor_2 % N
    # we'll run this until we get a pair of competitors that are actually different:
    while ((R[competitor_1[0], competitor_1[1]] == R[competitor_2t[0], competitor_2t[1]]) and
           (S[competitor_1[0], competitor_1[1]] == S[competitor_2t[0], competitor_2t[1]]) and
           (C[competitor_1[0], competitor_1[1]] == C[competitor_2t[0], competitor_2t[1]])):
        competitor_1, competitor_2 = competiroll(N)
        competitor_2t = competitor_2 % N
        # time passes:
        step += 1
#    print "competitor_1, competitor_2"
#    print competitor_1, competitor_2
    # here we produce torusified versions of the boards.
    # for signallers, we take both S_rad and C_rad around our competitors because,
    # signallers affect receptive && cooperating cells which affect our competitors
    S_sub = S[sp.arange(- S_rad - C_rad, S_rad + C_rad + 1)%N, :][:, sp.arange(- S_rad - C_rad, S_rad + C_rad + 1)%N]
    R_sub = R[sp.arange(- C_rad, C_rad + 1)%N, :][:, sp.arange(- C_rad, C_rad + 1)%N]
    C_sub = C[sp.arange(- C_rad, C_rad + 1)%N, :][:, sp.arange(- C_rad, C_rad + 1)%N]
#    print "S_sub.shape, R_sub.shape, C_sub.shape"
#    print S_sub.shape, R_sub.shape, C_sub.shape
    # we count how many signallers are within each cell's neighbourhood
    S_conv = sp.signal.convolve2d(S_sub, S_counter, mode='valid')
    # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
    cooping_cells = (C_sub == R_sub) == (S_conv > S_th)
    # how many cooperators around each competitor?
#    print "cooping_cells"
#    print cooping_cells.shape
#    print cooping_cells
    C_conv = sp.signal.convolve2d(cooping_cells, C_counter, mode='valid')
    # Public goods effect.
    # G for Goods
    G = (C_conv > C_th)
#    print "G.shape", G.shape
    # all cells for which the effect of goods is above threshold is True in G.
    # M for Metabolism
    cost_board = S_cost * S + R_cost * R + C_cost * C
    M = G * (1 - benefit) * cost_board
    # all false in G don't benefit from public goods (G^True flips values)
    M += (G^True) *  cost_board
    if M[competitor_1[0], competitor_1[1]] > M[competitor_2t[0], competitor_2t[1]]:
        C[competitor_1[0], competitor_1[1]] = C[competitor_2t[0], competitor_2t[1]]
        S[competitor_1[0], competitor_1[1]] = S[competitor_2t[0], competitor_2t[1]]
        R[competitor_1[0], competitor_1[1]] = R[competitor_2t[0], competitor_2t[1]]
    elif M[competitor_1[0], competitor_1[1]] == M[competitor_2t[0], competitor_2t[1]]:
        nothing #print 'buga'
    else:
        C[competitor_2t[0], competitor_2t[1]] = C[competitor_1[0], competitor_1[1]]
        S[competitor_2t[0], competitor_2t[1]] = S[competitor_1[0], competitor_1[1]]
        R[competitor_2t[0], competitor_2t[1]] = R[competitor_1[0], competitor_1[1]]    

def mutate():
    global C, R, S
    if sp.random.random()>0.1:
        coords = sp.random.randint(N, size=2)
        B = [C, R, S][sp.random.randint(3)]
        B[coords[0], coords[1]] = sp.random.randint(2)

def diffuse():
    global R, S, C
    m, n = sp.random.randint(N, size=2)
    m1, n1 = (m+1)%N, (n+1)%N
    if sp.random.randint(2):
        # Truth for clockwise
        for board in [R, S, C]:
            board[[m, m, m1, m1], [n, n1, n1, n]] = board[[m1, m, m, m1], [n, n, n1, n1]]
    else:
        for board in [R, S, C]:
            board[[m, m, m1, m1], [n, n1, n1, n]] = board[[m, m1, m1, m], [n1, n1, n, n]]

def sample_frequency():
    global S, R, C
    res = zeros(8, dtype='uint32')
    for i in sp.arange(8):
        res[i] = sum(bool8((not 1 & i) ^ S) +
                     bool8((not 2 & i) ^ R) +
                     bool8((not 4 & i) ^ C))
    return res

def sample_nhood():
    global S, R, C
    res = zeros((8, 8), dtype='uint32')
    
def sample():
    global sample_count, samples_frequency, samples_nhood
    if not step % steps_per_sample:
        samples_frequency[sample_count] = sample_count
        samples_nhood[sample_count] = sample_count
        sample_count += 1
            
def mainstuff():
    global while_count, step, data
    #print "while_count", while_count
    competition()
    diffuse()
    sample()
    step += 1
    while_count += 1
    print step

## process data

def imagify_data():
    global image, S, R, C
    ## package data
    data = sp.ndarray(shape=[N, N, 3])
    data[:, :, 0] = S
    data[:, :, 1] = R
    data[:, :, 2] = C
    resized_data = (255*data).repeat(4, axis=0).repeat(4, axis=1)
    image = pygame.surfarray.make_surface(resized_data)

## update display

def update_display():
    global image
    screen.blit(image, (0, 0))
    pygame.display.flip()

#clock = pygame.time.Clock()

# infinite loop

while_count = 0

while True:
    mainstuff() # the main algorithm
#    imagify_data()
#    update_display()
