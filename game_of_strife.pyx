#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

#import pyximport
#pyximport.install()
#import cython
import h5py
#import os
import time
import scipy as sp
import scipy.signal
import numpy as np
cimport numpy as np
#import pygame
#import pylab as pl
#import timeit
#import sys

cdef:
    # Board size
    int N = 10

    # we'll increase this by one every time two cells compete.
    np.uint32_t step_count = 0

    ## time keeping
    # number of generations the simulation will run
    # each generation is defined as the average number of steps for which
    # each cell on the board was in a competition once since last
    # generation (?)
    
    np.uint16_t generations = 10
    np.uint16_t steps_final = generations * N**2
    
    # Cost of gene expression
    np.double_t S_cost = 3
    np.double_t R_cost = 1
    np.double_t C_cost = 30
    np.double_t B_cost = 100  # B for Baseline

    # cooperation benefit, in ratio
    np.double_t benefit = 0.3
    
    # number of genotypes possible
    np.uint8_t genotype_num = 8
    
    # mutation per competition
    double mutation_rate = 1e-1

    # radius
    int S_rad = 1
    int C_rad = 1
    

    # diameter of the convolution matrix
    int S_len = 2 * S_rad + 1
    int C_len = 2 * C_rad + 1

    # the convolution matrix used to count neighbours
    np.ndarray S_kernel = sp.ones((S_len, S_len), np.uint16)
    np.ndarray C_kernel = sp.ones((C_len, C_len), np.uint16)
    
    ## neighbours effects' thresholds
    int S_th = 6
    # quorum threshold
    int C_th = 3
    # Cooperation threshold. Above it, public goods makes a difference.

    ## data sampling
    # we will take a frequency sample some number of times per generation
    np.uint32_t steps_per_gen = N ** 2
    np.uint8_t samples_per_gen = 1
    np.uint32_t samples_num = samples_per_gen * generations
    np.uint32_t sample_count = 0
    np.uint32_t steps_per_sample = sp.floor(1.0 * steps_per_gen / samples_per_gen)
    
    # We want to know the frequency of each genotype per generation
    np.ndarray samples_frequency = sp.empty((samples_num, genotype_num), dtype='uint32')
#    np.ndarray samples_nhood = sp.empty((samples_num, genotype_num, genotype_num))

#labels = ['Ignorant (crs)', 'Liar (crS)', 'Voyeur (cRs)', 'Lame (cRS)',
#          'Blunt (Crs)', 'Vain (CrS)', 'Shy (CRs)', 'Honest (CRS)']

    ## settings

#        cdef public np.uint_t N
#        cdef np.uint_t cell_num
#        cdef np.uint_t step_count
#        cdef np.uint_t generations

# pygame initialization

#pygame.init()
#screen = pygame.display.set_mode((N*4, N*4))
#pygame.display.set_caption("lets see")

## functions

#    @profile
cpdef int competition(np.ndarray[np.uint8_t, ndim=3] B) except -1:
    cdef:
        np.ndarray[np.uint16_t, ndim=1] c_pos_1, c_pos_2, c_pos_2t
        np.ndarray[np.uint16_t, ndim=1] s_r_range, s_c_range
        np.ndarray[np.uint16_t, ndim=1] rc_r_range, rc_c_range
        np.ndarray S_sub, R_sub, C_sub
        np.ndarray[np.uint8_t, ndim=2] S_conv, C_conv, cooping_cells
        np.ndarray G, M
        np.ndarray[np.double_t, ndim=2] S_cost_board, R_cost_board, C_cost_board
        np.ndarray[np.double_t, ndim=2] Total_cost_board
        np.uint16_t rl, rh, cl, ch
        np.double_t p0, p1
        np.ndarray[np.uint16_t, ndim=2] NEIGHBOUR_REL_POS = np.array([[0, -1], [0, 1], [-1, 0], [1, 0]], dtype=np.uint16)

    
    ##
    # Draw two adjacent positions.
    # We'll use relative positions to compute exact positions of 2nd competitor cell
    #cdef np.ndarray B = B
    c_pos_1 = sp.int16(sp.random.randint(N, size=2))
    c_pos_2 = c_pos_1 + NEIGHBOUR_REL_POS[sp.random.randint(4)]
    # c_pos_2's coordinates in a torus:
    c_pos_2t = c_pos_2 % N
    # two identical cells competing will result in two identical cells,
    # so we will return now with no further calculation of this competition.
    if (B[:, c_pos_1[0], c_pos_1[1]] == B[:, c_pos_2t[0], c_pos_2t[1]]).all():
        mutate(B, c_pos_1)
        return 1

    ## We will optimize by taking a sub array from each genotype array around the competitors.

    # rl, ch - row low, col high
    rl, rh = sp.sort([c_pos_1[0], c_pos_2[0]])
    cl, ch = sp.sort([c_pos_1[1], c_pos_2[1]])

    # For signallers, we take both S_rad and C_rad around our competitors because
    # signallers affect CR[Ss] cells which, with their public goods, affect our competitors
    s_r_range = sp.arange(rl - S_rad - C_rad, rh + S_rad + C_rad + 1, dtype=np.uint16) % N
    s_c_range = sp.arange(cl - S_rad - C_rad, ch + S_rad + C_rad + 1, dtype=np.uint16) % N
    rc_r_range = sp.arange(rl - C_rad, rh + C_rad + 1, dtype=np.uint16) % N
    rc_c_range = sp.arange(cl - C_rad, ch + C_rad + 1, dtype=np.uint16) % N

    S_sub = B[0, s_r_range, :][:, s_c_range]
    R_sub = B[1, rc_r_range, :][:, rc_c_range]
    C_sub = B[2, rc_r_range, :][:, rc_c_range]

#        raw_input(S_sub)
#        raw_input(S_sub.shape)
#        raw_input(R_sub)
#        raw_input(R_sub.shape)
#        raw_input(C_sub)
#        raw_input(C_sub.shape)
    #    print "S_sub.shape, R_sub.shape, C_sub.shape"
    #    print S_sub.shape, R_sub.shape, C_sub.shape

    # we count how many signallers are within each cell's neighbourhood
    #print S_sub.shape
    S_conv = np.uint8(sp.signal.convolve2d(S_sub, S_kernel, mode='valid'))

    # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
    # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
    cooping_cells = ((C_sub & R_sub) & (S_conv > S_th)) | (C_sub & (R_sub ^ True))
    # how many cooperators around each competitor?
    #    print "cooping_cells"
    #    print cooping_cells.shape
    #    print cooping_cells
    C_conv = sp.array(sp.signal.convolve2d(cooping_cells, C_kernel, mode='valid'), dtype=np.uint8)
    # Public goods effect.
    # G for Goods
    G = (C_conv > C_th)
    #    print "G.shape", G.shape
    # all cells for which the effect of goods is above threshold is True in G.
    # M for Metabolism
    S_cost_board = S_cost * B[0, sp.arange(rl, rh + 1) % N, :][:, sp.arange(cl, ch + 1) % N]
    R_cost_board = R_cost * B[1, sp.arange(rl, rh + 1) % N, :][:, sp.arange(cl, ch + 1) % N]
    C_cost_board = C_cost * B[2, sp.arange(rl, rh + 1) % N, :][:, sp.arange(cl, ch + 1) % N]
    Total_cost_board = S_cost_board + R_cost_board + C_cost_board + B_cost
    M = G * (1 - benefit) * Total_cost_board
    # all false in G don't benefit from public goods (G^True flips values)
    M += (G^True) *  Total_cost_board
    M = B_cost / M
    p0 = sp.rand() * M.item(0)
    p1 = sp.rand() * M.item(1)
    if c_pos_1[0] != c_pos_2[0]:
        if c_pos_1[0] > c_pos_2[0]:
            # their position is like this:
            # 2
            # 1
            if p0 > p1:
                #print "comp 2 wins"
                # competitor 2 wins
                endgame(B, c_pos_2t, c_pos_1)
            else:
                #print "comp 1 wins"
                # competitor 1 wins
                endgame(B, c_pos_1, c_pos_2t)
        else:
            # their position is like this:
            # 1
            # 2
            if p0 > p1:
                #print "comp 2 wins"
                # competitor 2 wins
                endgame(B, c_pos_2t, c_pos_1)
            else:
                #print "comp 1 wins"
                # competitor 1 wins
                endgame(B, c_pos_1, c_pos_2t)
    else:
        if p0 > p1:
            # their position is like this:
            # 2 1
            if M[0,0] > M[0,1]:
                #print "comp 2 wins"
                # competitor 2 wins
                endgame(B, c_pos_2t, c_pos_1)
            else:
                # competitor 1 wins
                endgame(B, c_pos_1, c_pos_2t)
        else:
            # their position is like this:
            # 1 2
            if p0 > p1:
                #print "comp 2 wins"
                # competitor 2 wins
                endgame(B, c_pos_2t, c_pos_1)
            else:
                #print "comp 1 wins"
                # competitor 1 wins
                endgame(B, c_pos_1, c_pos_2t)

cpdef endgame(np.ndarray[np.uint8_t, ndim = 3] B, np.ndarray[np.uint16_t, ndim = 1] winner, np.ndarray[np.uint16_t, ndim = 1] loser):
    copycell(B, winner, loser)
    mutate(B, loser)
    
cpdef copycell(np.ndarray[np.uint8_t, ndim = 3] B, np.ndarray[np.uint16_t, ndim = 1] orig, np.ndarray[np.uint16_t, ndim = 1] copy):
    B[:, copy[0], copy[1]] = B[:, orig[0], orig[1]]

cpdef mutate(np.ndarray[np.uint8_t, ndim = 3] B, np.ndarray[np.uint16_t, ndim = 1] pos):
    cdef np.uint8_t loci
    if sp.rand() < mutation_rate:
        loci = sp.random.randint(3)
        B[loci, pos[0], pos[1]] = not B[loci, pos[0], pos[1]]
#            B[loci, pos[0], pos[1]] = B[loci, pos[0], pos[1]] ^ B[loci, pos[0], pos[1]]

cpdef int diffuse(B) except -1:
    cdef np.uint16_t m, n, m1, n1
    m, n = sp.random.randint(N, size=2)
    m1, n1 = (m + 1) % N, (n + 1) % N
    if sp.random.rand() < 0.5:
        B[:, (m, m, m1, m1), (n, n1, n1, n)] = B[:, (m1, m, m, m1), (n, n, n1, n1)]
    else:
        B[:, (m, m, m1, m1), (n, n1, n1, n)] = B[:, (m1, m, m, m1), (n, n, n1, n1)]

cpdef int save(f) except -1:
    pass

#    def diffuse(self):
#        m, n = sp.random.randint(N, size=2)
#        m1, n1 = (m + 1) % N, (n + 1) % N
#        if sp.random.rand()<0.5:
#            # Truth for clockwise
#            for board in [R, S, C]:
#                board[(m, m, m1, m1), (n, n1, n1, n)] = board[(m1, m, m, m1), (n, n, n1, n1)]
#        else:
#            for board in [R, S, C]:
#                board[(m, m, m1, m1), (n, n1, n1, n)] = board[(m, m1, m1, m), (n1, n1, n, n)]

cpdef int sample(np.ndarray[np.uint8_t, ndim=3] B) except -1:
    global samples_frequency, sample_count
    cdef:
        np.uint8_t genotype
        np.uint32_t genotype_frequency
        
    for genotype in xrange(8):
        genotype_board = B[0] + 2 * B[1] + 4 * B[2] == genotype
        genotype_frequency = sp.sum(genotype_board)
#            for j in range(8):
#                genotype_board
#                samples_nhood[sample_count, genotype]
        samples_frequency[sample_count, genotype] = genotype_frequency
    sample_count += 1

cpdef int nextstep(np.ndarray[np.uint8_t, ndim=3] B) except -1:
    global step_count, steps_per_sample

    competition(B)
#        diffuse()
    if not step_count % steps_per_sample: sample(B)
    #print step_count
    step_count += 1

## process data

#cdef int stratificatied() except -1:
#    cdef:
#        np.ndarray[np.uint32_t, ndim=2] res
#    res = sp.empty((samples_num, genotype_num))
#    for i in xrange(genotype_num):
#        res[:,i] = sp.array([samples_frequency[:,i] + sp.sum(samples_frequency[:,:i])])
#    return res

#    def display_frequency_timeseries(self):
#        for i in range(8):
#            pl.plot(sp.arange(samples_num), samples_frequency[:,i], label=str(i), fillstyle='bottom')
cpdef int print_B(np.ndarray[np.uint8_t, ndim=3] B) except -1:
    print B


#def make_surface(image_data):
#    resized_data = (255 * image_data).repeat(4, axis=0).repeat(4, axis=1)
#    return pygame.surfarray.make_surface(resized_data)

### update display
#
#def update_display():
#    image
#    screen.blit(image, (0, 0))
#    pygame.display.flip()

#clock = pygame.time.Clock()

def go():
    # A cell can be Signalling and/or Receptive and/or Cooperative
    cdef:
        np.ndarray[np.uint8_t, ndim=3] B = sp.ones((3, N, N), dtype=np.uint8)
        np.double_t t
        np.uint32_t steps_0
#    raw_input('start')
    t = time.time()
    every = 10
    print t, step_count
    steps_0 = step_count
    while step_count < steps_final:
        nextstep(B)
        delta_t = time.time() - t
        if not step_count % 5000: print step_count
        if delta_t > every:
            savestrife(B)
            t = time.time()
            steps_delta = step_count - steps_0
            steps_0 = step_count
            print t, step_count, steps_delta
    print_B(B)
    print (time.time() - t)/generations/N/N*300*300*10000/60/60/24
def loadstrife(np.ndarray[np.uint8_t, ndim=3] B, filename='strife_in_a_jar.npz'):
    f = open(filename, 'rb')
    ff = np.load(f)
    N = ff['N']
    step_count = ff['step_count']
    generations = ff['generations']
    steps_final = ff['steps_final']
    genotype_num = ff['genotype_num']
    S_cost = ff['S_cost']
    R_cost = ff['R_cost']
    C_cost = ff['C_cost']
    B_cost = ff['B_cost']
    benefit = ff['benefit']
    mutation_rate = ff['mutation_rate']
    S_rad = ff['S_rad']
    C_rad = ff['C_rad']
    S_len = ff['S_len']
    C_len = ff['C_len']
    S_kernel = ff['S_kernel']
    C_kernel = ff['C_kernel']
    S_th = ff['S_th']
    C_th = ff['C_th']
    B = ff['B']
    steps_per_gen = ff['steps_per_gen']
    samples_per_gen = ff['samples_per_gen']
    samples_num = ff['samples_num']
    steps_per_sample = ff['steps_per_sample']
    sample_count = ff['sample_count']
    samples_frequency = ff['samples_frequency']
    return ff
#    samples_nhood = ff['samples_nhood']
def savestrife(np.ndarray[np.uint8_t, ndim=3] B, f='strife_in_a_jar.npz'):
    f = open(f, 'wb')
    ff = {}
    ff['N'] = N
    ff['step_count'] = step_count
    ff['generations'] = generations
    ff['steps_final'] = steps_final
    ff['genotype_num'] = genotype_num
    ff['S_cost'] = S_cost
    ff['R_cost'] = R_cost
    ff['C_cost'] = C_cost
    ff['B_cost'] = B_cost
    ff['benefit'] = benefit
    ff['mutation_rate'] = mutation_rate
    ff['S_rad'] = S_rad
    ff['C_rad'] = C_rad
    ff['S_len'] = S_len
    ff['C_len'] = C_len
    ff['S_kernel'] = S_kernel
    ff['C_kernel'] = C_kernel
    ff['S_th'] = S_th
    ff['C_th'] = C_th
    ff['B'] = B
    ff['steps_per_gen'] = steps_per_gen
    ff['samples_per_gen'] = samples_per_gen
    ff['samples_num'] = samples_num
    ff['steps_per_sample'] = steps_per_sample
    ff['sample_count'] = sample_count
    ff['samples_frequency'] = samples_frequency
#    ff['samples_nhood'] = samples_nhood
    np.savez(f, ff)
