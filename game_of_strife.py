#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import h5py
import os
import time
import scipy as sp
import scipy.signal
#import pygame
import pylab as pl
#import timeit
import sys

labels = ['Ignorant (csr)', 'Voyeur (csR)', 'Liar (cSr)', 'Lame (cSR)',
          'Blunt (Csr)', 'Shy (CsR)', 'Vain (CSr)', 'Honest (CSR)']

class game_of_strife:
    def __init__(self, N=10, generations=1):
        
        d = {}
        
        #########
        # model parameters
        #########
        
        #######
        # Cost of gene expression
        d['S_cost'] = sp.array(3)
        d['R_cost'] = sp.array(1)
        d['C_cost'] = sp.array(30)
        d['B_cost'] = sp.array(100)  # B for Baseline, basal, "basic metabolic burden"

        #####
        # benefit from cooperation. "reward factor" in the article.
        d['benefit'] = sp.array(0.9)

        ######
        # mutation per generation
        d['mutation_rate_r'] = sp.array(1e-4)
        d['mutation_rate_s'] = sp.array(1e-4)
        d['mutation_rate_c'] = sp.array(1e-4)

        ## neighbours effects' thresholds
        d['S_th'] = sp.array(6)
        # quorum threshold
        d['C_th'] = sp.array(3)
        # Cooperation threshold. Above it, public goods makes a difference.
        
        ## Probability of each single diffusion operation
        d['D'] = 0.2
        
        #####
        # settings
        #####
        
        d['NEIGHBOUR_REL_POS'] = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        # Board size
        d['N'] = sp.array(N)
        d['num_cells'] = sp.array(N ** 2)

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        #######
        # we'll increase step_count by one every time two cells compete.
        d['step_count'] = sp.array(0)


        d['generations'] = sp.array(generations)
        d['steps_final'] = sp.array(generations * d['num_cells'])

        #######
        # radius of Signal or Cooperation effects.
        d['S_rad'] = sp.array(1)
        d['C_rad'] = sp.array(1)

        # diameter of the convolution matrix
        diameter = lambda x: sp.array(2 * x + 1)
        d['S_len'] = diameter(d['S_rad'])
        d['C_len'] = diameter(d['C_rad'])

        # the convolution matrix used to count neighbours
        d['S_kernel'] = sp.ones((d['S_len'], d['S_len']))
        d['C_kernel'] = sp.ones((d['C_len'], d['C_len']))

        # A cell can be Signalling and/or Receptive and/or Cooperative
        R = sp.zeros((N, N), dtype='bool')
        S = sp.zeros((N, N), dtype='bool')
        d['B'] = sp.array([R, S])
        
        d['genotype_num'] = sp.array(4)
        
        ## data sampling
        # we will take a frequency sample some number of times per generation
        d['steps_per_gen'] = sp.array(N ** 2)
        d['samples_per_gen'] = sp.array(1)
        d['samples_num'] = sp.array(d['samples_per_gen'] * d['generations'])
        d['steps_per_sample'] = sp.array(sp.floor(1.0 * d['steps_per_gen'] / d['samples_per_gen']), dtype=sp.uint32)
        d['sample_count'] = sp.array(0)
        # We want to know the frequency of each genotype per generation
        d['samples_frequency'] = sp.empty((d['samples_num'], d['genotype_num']), dtype='int32')
        d['samples_nhood'] = sp.empty((d['samples_num'], d['genotype_num'], d['genotype_num']), dtype=sp.int32)
        
        self.__init_unpack_parameters__(d)
    
    def __init_unpack_parameters__(self, d):
        self.parameters = set()
        for key, val in d.items():
            setattr(self, key, val)
            self.parameters.add(key)
            
    ######
    ## functions
    ######    
    
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

        # we count how many signallers are within each cell's neighbourhood
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')

        # a cell will produce common goods if the signal in its neighborhood, that is compatible with its receptor, is above threshold.

        # R1 & (S1 > Sth) or  R2 & (S2 > Sth)
        type_2_above_threshold = S_conv > self.S_th
        type_2_cooping = R_sub & type_2_above_threshold
        type_1_cooping =  (R_sub | type_2_above_threshold) ^ True
        cooping_cells = type_1_cooping | type_2_cooping
        
        
#        # how many cooperators around each competitor?
#        C_conv = sp.signal.convolve2d(cooping_cells, self.C_kernel, mode='valid')
        
        
        # Public goods effect.
        # G for Goods.
        # Which are the cells that enjoy the effect of public goods?
        G = (cooping_cells > self.C_th)
        
        
        # all cells for which the effect of goods is above threshold is True in G.
        # M for Metabolism
        # Rock Paper Scissors => C.all() == True
        twocellpos_r, twocellpos_c = sp.arange(rl, rh + 1) % self.N, sp.arange(cl, ch + 1) % self.N
        twocellrange_c = sp.arange(cl, ch + 1)
        R_cost_board = self.R_cost * self.B[0, twocellpos_r, twocellpos_c]
        S_cost_board = self.S_cost * self.B[1, twocellpos_r, twocellpos_c]
        Total_cost_board = S_cost_board + R_cost_board + self.C_cost + self.B_cost
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
        if sp.rand() < self.mutation_rate_r:
            self.B[0, pos[0], pos[1]] = self.B[0, pos[0], pos[1]] ^ True
        if sp.rand() < self.mutation_rate_s:
            self.B[1, pos[0], pos[1]] = self.B[1, pos[0], pos[1]] ^ True
        if sp.rand() < self.mutation_rate_c:
            self.B[2, pos[0], pos[1]] = self.B[2, pos[0], pos[1]] ^ True

    def diffuse(self):
        m, n = sp.random.randint(self.N, size=2)
        m1, n1 = (m + 1) % self.N, (n + 1) % self.N
        if sp.random.rand() < 0.5:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]
        else:
            self.B[:, (m, m, m1, m1), (n, n1, n1, n)] = self.B[:, (m1, m, m, m1), (n, n, n1, n1)]
        
    def save_h5(self, fname):
        with h5py.File(fname) as ff:
            for key in self.parameters:
                try:
                    print key, type(getattr(self, key))
                    ff[key] = getattr(self, key)
                except:
                    ff[key][...] = getattr(self, key)
    
    def load_h5(self, fname):
        with h5py.File(fname) as ff:
            d = { key : val[...] for key, val in ff.items() }
            self.__init_unpack_parameters__(d)
    
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
        for i in range(sp.floor(self.N ** 2 * self.D):
            self.diffuse()
        if not self.step_count % self.steps_per_sample: self.sample()
        print self.step_count
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
    every = 10
    print "t: %(t)f, steps thus far: %(steps)d" % {'t': t, 'steps': a.step_count}
    steps_0 = a.step_count
    while a.step_count < a.steps_final:
        print a.step_count
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            t = time.time()
            a.save_h5(fname)
            steps_delta = a.step_count - steps_0
            steps_0 = a.step_count
            eta = 1.0 * delta_t / steps_delta * (a.steps_final - a.step_count)
            print "t: %(t)f, approx. time to fin: %(eta)f" % {'t': t, 'eta': eta}
            print "steps taken = %(step_count)s, steps since last save = %(steps_delta)s" % {'step_count': a.step_count, 'steps_delta': steps_delta}
            sys.exit(1)

if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except IndexError as e:
        print "Usage: %(me)s [datafilename]" % { 'me' : sys.argv[0] }
        raise
    a = game_of_strife(N=10, generations=1)
    if os.path.exists(fname):
        a.load_h5(fname)
        go(a)
    else:
        a = game_of_strife(N=100, generations=10000)
        a.save_h5(fname)
        go(a)
    a.save_h5(fname)
    sys.exit(0)