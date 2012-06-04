#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import cython
import h5py
import copy
import os
import gc
import time
import pickle
import scipy as sp
import scipy.signal
#import pygame
import pylab as pl
#import timeit
#import sys

labels = ['Ignorant (csr)', 'Voyeur (csR)', 'Liar (cSr)', 'Lame (cSR)',
          'Blunt (Csr)', 'Shy (CsR)', 'Vain (CSr)', 'Honest (CSR)']

class game_of_strife:
    def __init__(self, N=10, generations=1):
        ## settings
        
        d = {}
        
        d['data_filename'] = 'data.npz'
        
        d['NEIGHBOUR_REL_POS'] = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        # Board size
        d['N'] = N
        d['num_cells'] = N ** 2

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        #######
        # we'll increase step_count by one every time two cells compete.
        d['step_count'] = 0


        d['generations'] = generations
        d['steps_final'] = generations * d['num_cells']

        #######
        # number of genotypes possible
        d['genotype_num'] = 8

        #######
        # Cost of gene expression
        d['S_cost'] = 3
        d['R_cost'] = 1
        d['C_cost'] = 30
        d['B_cost'] = 100  # B for Baseline

        #####
        # cooperation benefit, in ratio
        d['benefit'] = 0.3

        ######
        # mutation per generation
        d['mutation_rate'] = 1e-4

        #######
        # radius of Signal or Cooperation effects.
        d['S_rad'] = 1
        d['C_rad'] = 1

        # diameter of the convolution matrix
        diameter = lambda x: 2 * x + 1
        d['S_len'] = diameter(d['S_rad'])
        d['C_len'] = diameter(d['C_rad'])

        # the convolution matrix used to count neighbours
        d['S_kernel'] = sp.ones((d['S_len'], d['S_len']))
        d['C_kernel'] = sp.ones((d['C_len'], d['C_len']))

        ## neighbours effects' thresholds
        d['S_th'] = 6
        # quorum threshold
        d['C_th'] = 3
        # Cooperation threshold. Above it, public goods makes a difference.

        # A cell can be Signalling and/or Receptive and/or Cooperative
        R = sp.zeros((N, N), dtype='bool')
        S = sp.zeros((N, N), dtype='bool')
        C = sp.ones((N, N), dtype='bool')
        d['B'] = sp.array([R, S, C])

        ## data sampling
        # we will take a frequency sample some number of times per generation
        d['steps_per_gen'] = N ** 2
        d['samples_per_gen'] = 1
        d['samples_num'] = d['samples_per_gen'] * d['generations']
        d['steps_per_sample'] = sp.uint32(sp.floor(1.0 * d['steps_per_gen'] / d['samples_per_gen']))
        d['sample_count'] = 0
        # We want to know the frequency of each genotype per generation
        d['samples_frequency'] = sp.empty((d['samples_num'], d['genotype_num']), dtype='int32')
        d['samples_nhood'] = sp.empty((d['samples_num'], d['genotype_num'], d['genotype_num']))
        
        self.__init_unpack_parameters__(d)
    
    def __init_unpack_parameters__(self, d):
        self.parameters = []
        for key, val in d.items():
            setattr(self, key, val)
            self.parameters.append(key)
            

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
#        twosort = lambda x, y: (x, y) if x < y else (y, x)
        rl, rh = (c_pos_1[0], c_pos_2[0]) if c_pos_1[0] < c_pos_2[0] else (c_pos_2[0], c_pos_1[0])
        cl, ch = (c_pos_1[1], c_pos_2[1]) if c_pos_1[1] < c_pos_2[1] else (c_pos_2[1], c_pos_1[1])
        
        # For signallers, we take both S_rad and C_rad around our competitors because
        # signallers affect CR[Ss] cells which, with their public goods, affect our competitors
        s_row_range = sp.arange(rl - self.S_rad - self.C_rad, rh + self.S_rad + self.C_rad + 1) % self.N
        s_col_range = sp.arange(cl - self.S_rad - self.C_rad, ch + self.S_rad + self.C_rad + 1) % self.N
        rc_row_range = sp.arange(rl - self.C_rad, rh + self.C_rad + 1) % self.N
        rc_col_range = sp.arange(cl - self.C_rad, ch + self.C_rad + 1) % self.N

        R_sub = self.B[0, rc_row_range, :][:, rc_col_range]
        S_sub = self.B[1, s_row_range, :][:, s_col_range]
        C_sub = sp.ones((2 * self.C_rad + sp.absolute(rl-rh), 2 * self.C_rad + sp.absolute(cl-ch)), dtype='bool')

        # we count how many signallers are within each cell's neighbourhood
        #print S_sub.shape
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')

        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
        # We'll use the C.all() == True version for the Rock Paper Sciccors version.

        # TODO: Reactivate when switching back to the general simulation.        
        # commenting the more general version of an unfixed C allele.
#        cooping_cells = ((C_sub & R_sub) & (S_conv > self.S_th)) | (C_sub & (R_sub ^ True))

        # this version assumes C is always the True allele.
        cooping_cells = (R_sub & (S_conv > self.S_th)) | (R_sub ^ True)
        
        
        # how many cooperators around each competitor?
        C_conv = sp.signal.convolve2d(cooping_cells, self.C_kernel, mode='valid')
        
        
        # Public goods effect.
        # G for Goods.
        # Which are the cells that enjoy the effect of public goods?
        G = (C_conv > self.C_th)
        
        
        # all cells for which the effect of goods is above threshold is True in G.
        # M for Metabolism
        # Rock Paper Scissors => C.all() == True
        twocellpos_r, twocellpos_c = sp.arange(rl, rh + 1) % self.N, sp.arange(cl, ch + 1) % self.N
        twocellrange_c = sp.arange(cl, ch + 1)
        R_cost_board = self.R_cost * self.B[0, twocellpos_r, twocellpos_c]
        S_cost_board = self.S_cost * self.B[1, twocellpos_r, twocellpos_c]
        C_cost_board = self.C_cost * sp.array([True, True], dtype='bool')
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
        self.copycell(winner, loser)
        self.mutate(loser)
        
    def copycell(self, orig, copy):
        self.B[:, copy[0], copy[1]] = self.B[:, orig[0], orig[1]]

    def mutate(self, pos):
        if sp.rand() < self.mutation_rate:
            loci = sp.random.randint(2) # We will only mutate S and R, never C. C stays uppercase. R = 0, S = 1, C = 2.
            self.B[loci, pos[0], pos[1]] = self.B[loci, pos[0], pos[1]] ^ True

    def diffuse(self):
        m, n = sp.random.randint(self.N, size=2)
        m1, n1 = (m + 1) % self.N, (n + 1) % self.N
        if sp.random.rand() < 0.5:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]
        else:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]
        
    def savestrife(self, fname=None):
        fname = fname if fname else self.data_filename
        ff = {}
        for key in self.parameters:
            ff[key] = getattr(self, key)
        sp.savez(fname, ff)
    
    def loadstrife(self, fname=None):
        fname = fname if fname else self.data_filename
        ff = sp.load(fname)['arr_0'].tolist()
        self.__init_unpack_parameters__(ff)
    
    def sample(self):
        joint_board = self.B[0] + 2 * self.B[1] + 4 * self.B[2]
        for genotype in range(8):
            genotype_board = joint_board == genotype
            genotype_frequency = sp.sum(genotype_board)
            # neighbours_genotype
            for nh_genotype in range(8):
                nh_board = joint_board == nh_genotype
                nh_genotype_count = sp.signal.convolve2d(nh_board, sp.ones((3,3)), mode='same', boundary='wrap')
                nh_genotype_count_of_genotype = sp.sum(nh_genotype_count * genotype_board, dtype=sp.int32)
                self.samples_nhood[self.sample_count, genotype, nh_genotype] = nh_genotype_count_of_genotype
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

def go(a):
    t = time.time()
    every = 1
    print t, a.step_count, "yo"
    steps_0 = a.step_count
    while a.step_count < a.steps_final:
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            t = time.time()
            a.savestrife()
            steps_delta = a.step_count - steps_0
            steps_0 = a.step_count
            print t, 1.0 * delta_t / steps_delta * (a.steps_final - a.step_count)
            print t, a.step_count, steps_delta

if __name__ == '__main__':
    a = game_of_strife(N=10, generations=1)
    if os.path.exists(a.data_filename):
        a.loadstrife()
        go(a)
    else:
        a = game_of_strife(N=300, generations=1000)
        go(a)
    a.savestrife()
    
    print a.B
    print a.samples_nhood
    pl.hold(True)
    for genotype in range(8):
        pl.plot(a.samples_frequency[:, genotype], label=labels[genotype])
    pl.show()
    for genotype in range(8):
        pl.plot(a.samples_nhood[:, genotype, :], label=labels[genotype])
    pl.show()