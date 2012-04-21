#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import scipy as sp
import scipy.signal
import pygame

class game_of_strife:
    def __init__(self, N=10):
        ## settings
        
        # Board size
        self.N = N
        self.cell_num = self.N**2
        
        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which each cell
        # on the board was in a competition once since last generation (?)

        # we'll increase this by one every time two cells compete.
        self.step = 0

        self.generations = 5
        self.steps_final = self.generations * self.cell_num

        # number of genotypes possible
        self.genotype_num = 8

        # Cost of gene expression

        self.S_cost = 3
        self.R_cost = 8
        self.C_cost = 30
        self.B_cost = 100 # B for Baseline

        # cooperation benefit, in ratio
        self.benefit = 0.3
        
        # mutation per generation
        self.mutationrate = 0.001
        
        # radius
        self.S_rad = 1
        self.C_rad = 1

        # diameter of the convolution matrix
        diameter = lambda x: 2 * x + 1
        self.S_len = diameter(self.S_rad)
        self.C_len = diameter(self.C_rad)

        # the convolution matrix used to count neighbours
        self.S_kernel = sp.ones((self.S_len, self.S_len))
        self.C_kernel = sp.ones((self.C_len, self.C_len)) # convolution matrix used to count cells that produce public goods

        # neighbours effects' thresholds
        self.S_th = 3 # quorum threshold
        self.C_th = 3 # Cooperation threshold. Above it, public goods makes a difference.

        # A cell can be Signalling and/or Receptive and/or Cooperative
        self.S = sp.rand(self.N, self.N) < 0.5
        self.R = sp.rand(self.N, self.N) < 0.5
        self.C = sp.rand(self.N, self.N) < 0.5

        ## data sampling
        # we will take a frequency sample some number of times per generation
        self.steps_per_gen = self.N**2
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

    def competiroll(self):
        #"""draw two competitor positions stuff"""
        # We'll use relative positions to compute exact positions of 2nd competitor cell
        NEIGHBOUR_ROW = sp.array([-1,  0,  1, -1,  0,  1, -1,  1])
        NEIGHBOUR_COL = sp.array([-1, -1, -1,  1,  1,  1,  0,  0])
        NEIGHBOUR_REL_POS = sp.array(zip(NEIGHBOUR_ROW, NEIGHBOUR_COL))
        c1 = sp.random.randint(self.N, size=2)
        c2 = c1 + NEIGHBOUR_REL_POS[sp.random.randint(8, size=1)[0]]
        return c1, c2

    def competition(self):
        # compete
        competitor_1, competitor_2 = self.competiroll()
        # competitor_2's coordinates in a torus:
        competitor_2t = competitor_2 % self.N
        # two identical cells competing will result in two identical cells,
        # so no need for further calculation.
        if ((self.R[competitor_1[0], competitor_1[1]] == self.R[competitor_2t[0], competitor_2t[1]]) and
               (self.S[competitor_1[0], competitor_1[1]] == self.S[competitor_2t[0], competitor_2t[1]]) and
               (self.C[competitor_1[0], competitor_1[1]] == self.C[competitor_2t[0], competitor_2t[1]])):
            return
        # print "competitor_1, competitor_2"
        # print competitor_1, competitor_2
        # here we produce torusified versions of the boards.
        # for signallers, we take both S_rad and self.C_rad around our competitors because,
        # signallers affect receptive && cooperating cells which affect our competitors
        S_sub = self.S[sp.arange(- self.S_rad - self.C_rad, self.S_rad + self.C_rad + 1)%self.N, :][:, sp.arange(- self.S_rad - self.C_rad, self.S_rad + self.C_rad + 1)%self.N]
        R_sub = self.R[sp.arange(- self.C_rad, self.C_rad + 1)%self.N, :][:, sp.arange(- self.C_rad, self.C_rad + 1)%self.N]
        C_sub = self.C[sp.arange(- self.C_rad, self.C_rad + 1)%self.N, :][:, sp.arange(- self.C_rad, self.C_rad + 1)%self.N]
        #    print "S_sub.shape, R_sub.shape, C_sub.shape"
        #    print S_sub.shape, R_sub.shape, C_sub.shape
        # we count how many signallers are within each cell's neighbourhood
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')
        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        cooping_cells = (C_sub == R_sub) == (S_conv > self.S_th)
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
        cost_board = self.B_cost + self.S_cost * self.S + self.R_cost * self.R + self.C_cost * self.C
        M = G * (1 - self.benefit) * cost_board
        # all false in G don't benefit from public goods (G^True flips values)
        M += (G^True) *  cost_board
        if M[competitor_1[0], competitor_1[1]] > M[competitor_2t[0], competitor_2t[1]]:
            self.C[competitor_1[0], competitor_1[1]] = self.C[competitor_2t[0], competitor_2t[1]]
            self.S[competitor_1[0], competitor_1[1]] = self.S[competitor_2t[0], competitor_2t[1]]
            self.R[competitor_1[0], competitor_1[1]] = self.R[competitor_2t[0], competitor_2t[1]]
#        elif M[competitor_1[0], competitor_1[1]] == M[competitor_2t[0], competitor_2t[1]]:
#            nothing #print 'buga'
        else:
            self.C[competitor_2t[0], competitor_2t[1]] = self.C[competitor_1[0], competitor_1[1]]
            self.S[competitor_2t[0], competitor_2t[1]] = self.S[competitor_1[0], competitor_1[1]]
            self.R[competitor_2t[0], competitor_2t[1]] = self.R[competitor_1[0], competitor_1[1]]    

    def mutate(self):
        if sp.random.random() < 0.05:
            coords = sp.random.randint(self.N, size=2)
            B = [self.C, self.R, self.S][sp.random.randint(3)]
            B[coords[0], coords[1]] = sp.random.randint(2)

    def diffuse(self):
        m, n = sp.random.randint(self.N, size=2)
        m1, n1 = (m+1)%self.N, (n+1)%self.N
        if sp.random.randint(2):
            # Truth for clockwise
            for board in [self.R, self.S, self.C]:
                board[[m, m, m1, m1], [n, n1, n1, n]] = board[[m1, m, m, m1], [n, n, n1, n1]]
        else:
            for board in [self.R, self.S, self.C]:
                board[[m, m, m1, m1], [n, n1, n1, n]] = board[[m, m1, m1, m], [n1, n1, n, n]]

    def sample(self):
        for genotype in range(8):
            genotype_board = self.S + 2 * self.R + 4 * self.C == genotype
            genotype_frequency = sp.sum(genotype_board)
#            for j in range(8):
#                genotype_board
#                self.samples_nhood[self.sample_count, genotype]
            self.samples_frequency[self.sample_count, genotype] = genotype_frequency
        self.sample_count += 1
                
    def nextstep(self):
        self.competition()
        self.mutate()
        self.diffuse()
        if not self.step % self.steps_per_sample: self.sample()
        #print self.step
        self.step += 1

    ## process data
    
    def stratificatied(self):
        res = sp.empty((self.sample_num, self.genotype_num))
        for i in range(self.genotype_num):
            res[:,i] = sp.array([self.samples_frequency[:,i] + sp.sum(self.samples_frequency[:,:i])])
        return res
            
    def imagify_data(self):
        ## package boards' data into a displayable array.
        return sp.array([self.S, self.R, self.C])
    
    def display_frequency_timeseries(self):
        for i in range(8):
            sp.plot(sp.arange(self.samples_num), self.samples_frequency[:,i], label=str(i), fillstyle='bottom')

def make_surface(image_data):
        resized_data = (255 * image_data).repeat(4, axis=0).repeat(4, axis=1)
        return pygame.surfarray.make_surface(resized_data)

### update display
#
#def update_display():
#    self.image
#    self.screen.blit(image, (0, 0))
#    self.pygame.display.flip()

#clock = pygame.time.Clock()

if __name__ == '__main__':
    a = game_of_strife()
    for i in range(a.steps_final):
        a.nextstep()
    #    imagify_data()
    #    update_display()
    #    print while_count
    sp.plot(a.samples_num, a.samples_frequency)
    raw_input('Cool, eh?')