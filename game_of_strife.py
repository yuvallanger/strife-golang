#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import os
import gc
import time
import pickle
import scipy as sp
import scipy.signal
#import pygame
#import pylab as pl
#import timeit
#import sys

class game_of_strife:
    def __init__(self, N=10, generations=1):
        ## settings
        
        self.NEIGHBOUR_REL_POS = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        # Board size
        self.N = N
        self.cell_num = self.N ** 2

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        # we'll increase this by one every time two cells compete.
        self.step_count = 0

        self.generations = generations
        self.steps_final = self.generations * self.cell_num

        # number of genotypes possible
        self.genotype_num = 8

        # Cost of gene expression

        self.S_cost = 3
        self.R_cost = 1
        self.C_cost = 30
        self.B_cost = 100  # B for Baseline

        # cooperation benefit, in ratio
        self.benefit = 0.3

        # mutation per generation
        self.mutation_rate = 1e-4

        # radius
        self.S_rad = 1
        self.C_rad = 1

        # diameter of the convolution matrix
        diameter = lambda x: 2 * x + 1
        self.S_len = diameter(self.S_rad)
        self.C_len = diameter(self.C_rad)

        # the convolution matrix used to count neighbours
        self.S_kernel = sp.ones((self.S_len, self.S_len))
        self.C_kernel = sp.ones((self.C_len, self.C_len))

        ## neighbours effects' thresholds
        self.S_th = 6
        # quorum threshold
        self.C_th = 3
        # Cooperation threshold. Above it, public goods makes a difference.

        # A cell can be Signalling and/or Receptive and/or Cooperative
        self.B = sp.ones((3, self.N, self.N), dtype='bool')
#        self.S = self.B[0]
#        self.R = self.B[1]
#        self.C = self.B[2]

        ## data sampling
        # we will take a frequency sample some number of times per generation
        self.steps_per_gen = self.N ** 2
        self.samples_per_gen = 1
        self.samples_num = self.samples_per_gen * self.generations
        self.steps_per_sample = sp.uint32(sp.floor(1.0 * self.steps_per_gen / self.samples_per_gen))
        self.sample_count = 0
        # We want to know the frequency of each genotype per generation
        self.samples_frequency = sp.empty((self.samples_num, self.genotype_num), dtype='int32')
        self.samples_nhood = sp.empty((self.samples_num, self.genotype_num, self.genotype_num))

        # pygame initialization

        #pygame.init()
        #screen = pygame.display.set_mode((N*4, N*4))
        #pygame.display.set_caption("lets see")

    ## functions
    
#    @profile
    def competition(self):
        ##
        # Draw two adjacent positions.
        # We'll use relative positions to compute exact positions of 2nd competitor cell
        c_pos_1 = sp.random.randint(self.N, size=2)
        c_pos_2 = c_pos_1 + self.NEIGHBOUR_REL_POS[sp.random.randint(4)]
        # c_pos_2's coordinates in a torus:
        c_pos_2t = c_pos_2 % self.N
        # two identical cells competing will result in two identical cells,
        # so we will return now with no further calculation of this competition.
        if (self.B[:, c_pos_1[0], c_pos_1[1]] == self.B[:, c_pos_2t[0], c_pos_2t[1]]).all():
            self.mutate(c_pos_1)
            return

        ## We will optimize by taking a sub array from each genotype array around the competitors.

        # rl, ch - row low, col high
        rl, rh = sp.sort([c_pos_1[0], c_pos_2[0]])
        cl, ch = sp.sort([c_pos_1[1], c_pos_2[1]])

        # For signallers, we take both S_rad and C_rad around our competitors because
        # signallers affect CR[Ss] cells which, with their public goods, affect our competitors
        s_r_range = sp.arange(rl - self.S_rad - self.C_rad, rh + self.S_rad + self.C_rad + 1) % self.N
        s_c_range = sp.arange(cl - self.S_rad - self.C_rad, ch + self.S_rad + self.C_rad + 1) % self.N
        rc_r_range = sp.arange(rl - self.C_rad, rh + self.C_rad + 1) % self.N
        rc_c_range = sp.arange(cl - self.C_rad, ch + self.C_rad + 1) % self.N

        S_sub = self.B[0, s_r_range, :][:, s_c_range]
        R_sub = self.B[1, rc_r_range, :][:, rc_c_range]
        C_sub = self.B[2, rc_r_range, :][:, rc_c_range]

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
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')

        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
        cooping_cells = ((C_sub & R_sub) & (S_conv > self.S_th)) | (C_sub & (R_sub ^ True))
        # how many cooperators around each competitor?
        #    print "cooping_cells"
        #    print cooping_cells.shape
        #    print cooping_cells
        C_conv = sp.signal.convolve2d(cooping_cells, self.C_kernel, mode='valid')
        # Public goods effect.
        # G for Goods
        G = (C_conv > self.C_th)
        #    print "G.shape", G.shape
        # all cells for which the effect of goods is above threshold is True in G.
        # M for Metabolism
        S_cost_board = self.S_cost * self.B[0, sp.arange(rl, rh + 1) % self.N, :][:, sp.arange(cl, ch + 1) % self.N]
        R_cost_board = self.R_cost * self.B[1, sp.arange(rl, rh + 1) % self.N, :][:, sp.arange(cl, ch + 1) % self.N]
        C_cost_board = self.C_cost * self.B[2, sp.arange(rl, rh + 1) % self.N, :][:, sp.arange(cl, ch + 1) % self.N]
        Total_cost_board = S_cost_board + R_cost_board + C_cost_board + self.B_cost
        M = G * (1 - self.benefit) * Total_cost_board
        # all false in G don't benefit from public goods (G^True flips values)
        M += (G^True) *  Total_cost_board
        M = self.B_cost / M
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
                    self.endgame(c_pos_2t, c_pos_1)
                else:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    self.endgame(c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1
                # 2
                if p0 > p1:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    self.endgame(c_pos_2t, c_pos_1)
                else:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    self.endgame(c_pos_1, c_pos_2t)
        else:
            if p0 > p1:
                # their position is like this:
                # 2 1
                if M[0,0] > M[0,1]:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    self.endgame(c_pos_2t, c_pos_1)
                else:
                    # competitor 1 wins
                    self.endgame(c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1 2
                if p0 > p1:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    self.endgame(c_pos_2t, c_pos_1)
                else:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    self.endgame(c_pos_1, c_pos_2t)
    
    def endgame(self, winner, loser):
        self.copycell(copy=loser, orig=winner)
        self.mutate(loser)
        
    def copycell(self, orig, copy):
        self.B[:, copy[0], copy[1]] = self.B[:, orig[0], orig[1]]

    def mutate(self, pos):
        if sp.rand() < self.mutation_rate:
            loci = sp.random.randint(3)
            self.B[loci, pos[0], pos[1]] = self.B[loci, pos[0], pos[1]] ^ self.B[loci, pos[0], pos[1]]

    def diffuse(self):
        m, n = sp.random.randint(self.N, size=2)
        m1, n1 = (m + 1) % self.N, (n + 1) % self.N
        if sp.random.rand()<0.5:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]
        else:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]

#    def diffuse(self):
#        m, n = sp.random.randint(self.N, size=2)
#        m1, n1 = (m + 1) % self.N, (n + 1) % self.N
#        if sp.random.rand()<0.5:
#            # Truth for clockwise
#            for board in [self.R, self.S, self.C]:
#                board[(m, m, m1, m1), (n, n1, n1, n)] = board[(m1, m, m, m1), (n, n, n1, n1)]
#        else:
#            for board in [self.R, self.S, self.C]:
#                board[(m, m, m1, m1), (n, n1, n1, n)] = board[(m, m1, m1, m), (n1, n1, n, n)]

    def sample(self):
        for genotype in range(8):
            genotype_board = self.B[0] + 2 * self.B[1] + 4 * self.B[2] == genotype
            genotype_frequency = sp.sum(genotype_board)
#            for j in range(8):
#                genotype_board
#                self.samples_nhood[self.sample_count, genotype]
            self.samples_frequency[self.sample_count, genotype] = genotype_frequency
        self.sample_count += 1

    def nextstep(self):
        self.competition()
#        self.diffuse()
        if not self.step_count % self.steps_per_sample: self.sample()
        #print self.step_count
        self.step_count += 1

    ## process data

    def stratificatied(self):
        res = sp.empty((self.samples_num, self.genotype_num))
        for i in range(self.genotype_num):
            res[:,i] = sp.array([self.samples_frequency[:,i] + sp.sum(self.samples_frequency[:,:i])])
        return res
            
    def imagify_data(self):
        ## package boards' data into a displayable array.
        return sp.array([self.S, self.R, self.C])
    
#    def display_frequency_timeseries(self):
#        for i in range(8):
#            pl.plot(sp.arange(self.samples_num), self.samples_frequency[:,i], label=str(i), fillstyle='bottom')
    


#def make_surface(image_data):
#    resized_data = (255 * image_data).repeat(4, axis=0).repeat(4, axis=1)
#    return pygame.surfarray.make_surface(resized_data)

### update display
#
#def update_display():
#    self.image
#    self.screen.blit(image, (0, 0))
#    self.pygame.display.flip()

#clock = pygame.time.Clock()

def picklize(a):
    f = open('strife_in_a_jar.pk', 'wb')
    pickle.dump(a, f, pickle.HIGHEST_PROTOCOL)
    f.close()

def go(a):
    t = time.time()
    every = 1200
    print t, a.step_count
    while a.step_count < a.steps_final:
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            picklize(a)
            t = time.time()
            gc.collect()
            print t, a.step_count

if __name__ == '__main__':
    if os.path.exists(r'strife_in_a_jar.pk'):
        a = pickle.load(open('strife_in_a_jar.pk', 'rb'))
        go(a)
    else:
        a = game_of_strife(N=300, generations=10000)
        go(a)
    picklize(a)
        
    #    imagify_data()
    #    update_display()
    #    print while_count
#    t = sp.arange(a.samples_num)
#    print t, a.samples_frequency
#    pl.hold(True)
#    pl.plot(t, a.samples_frequency)
#    pl.savefig('test.png')
#    raw_input('Cool, eh?')
