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
import scipy.weave
#import pygame
import pylab as pl
#import timeit
import sys
import signal

labels = ['Ignorant (csr)', 'Voyeur (csR)', 'Liar (cSr)', 'Lame (cSR)',
          'Blunt (Csr)', 'Shy (CsR)', 'Vain (CSr)', 'Honest (CSR)']

class game_of_strife:
    def __init__(self, fname, N=10, generations=1):
        
        d = {}
        
        ########
        # filname
        d['fname'] = fname

        #########
        # model parameters
        #########
        
        #######
        # Cost of gene expression
        d['S_cost'] = sp.uint64(3)
        d['R_cost'] = sp.uint64(1)
        d['C_cost'] = sp.uint64(30)
        d['B_cost'] = sp.uint64(100)  # B for Baseline, basal, "basic metabolic burden"

        #####
        # benefit from cooperation. "reward factor" in the article.
        d['benefit'] = sp.float64(0.9)

        ######
        # mutation per generation
        d['mutation_rate_r'] = sp.float64(1e-4)
        d['mutation_rate_s'] = sp.float64(1e-4)
        d['mutation_rate_c'] = sp.float64(1e-4)

        ## neighbours effects' thresholds
        d['S_th'] = sp.uint64(3)
        # quorum threshold
        d['C_th'] = sp.uint64(3)
        # Cooperation threshold. Above it, public goods makes a difference.
        
        ## Probability of each single diffusion operation
        d['D'] = sp.float64(0.5)
        
        #####
        # settings
        #####
        
        d['NEIGHBOUR_REL_POS'] = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        # Board size
        d['N'] = N
        d['num_cells'] = sp.uint64(N ** 2)

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        #######
        # we'll increase step_count by one every time two cells compete.
        d['step_count'] = sp.uint64(0)


        d['generations'] = sp.uint64(generations)
        d['steps_final'] = sp.uint64(generations * d['num_cells'])

        #######
        # radius of Signal or Cooperation effects.
        d['S_rad'] = sp.uint64(1)
        d['C_rad'] = sp.uint64(1)

        # diameter of the convolution matrix
        diameter = lambda x: sp.uint64(2 * x + 1)
        d['S_len'] = diameter(d['S_rad'])
        d['C_len'] = diameter(d['C_rad'])

        # the convolution matrix used to count neighbours
        d['S_kernel'] = sp.ones((d['S_len'], d['S_len']))
        d['C_kernel'] = sp.ones((d['C_len'], d['C_len']))

        # A cell can be Signalling and/or Receptive and/or Cooperative
        R = sp.rand(N, N) > 0.5
        S = sp.rand(N, N) > 0.5
        C = sp.rand(N, N) > 0.5
        d['B'] = sp.array([R, S, C])
        
        d['genotype_num'] = sp.uint64(8)
        
        ## data sampling
        # we will take a frequency sample some number of times per generation
        d['steps_per_gen'] = sp.uint64(N ** 2)
        d['samples_per_gen'] = sp.uint64(1)
        d['samples_num'] = sp.uint64(d['samples_per_gen'] * d['generations'])
        d['samples_board_num'] =sp.uint64(d['samples_num'] / 10)
        d['steps_per_sample'] = sp.uint64(sp.floor(1.0 * d['steps_per_gen'] / d['samples_per_gen']), dtype=sp.uint64)
        d['steps_per_board_sample'] = sp.uint64(10 * d['steps_per_sample'])
        d['sample_count'] = sp.uint64(0)
        # We want to know the frequency of each genotype per generation
        d['samples_frequency'] = sp.empty((d['samples_num'], d['genotype_num']), dtype='int32')
        d['samples_nhood'] = sp.empty((d['samples_num'], d['genotype_num'], d['genotype_num']), dtype=sp.uint64)
        d['samples_board'] = sp.empty((d['samples_board_num'], 3, N, N), dtype=sp.uint64)
        
        self.__init_unpack_parameters__(d)
    
    def __init_unpack_parameters__(self, d):
        self.parameters = set()
        for key, val in d.items():
            setattr(self, key, val)
            self.parameters.add(key)
            
    ######
    ## functions
    ######    
    
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
        C_sub = self.B[2, rc_row_range, :][:, rc_col_range]

        # we count how many signallers are within each cell's neighbourhood
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')

        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
        cooping_cells = ((C_sub & R_sub) & (S_conv >= self.S_th)) | (C_sub & (R_sub ^ True))
        
        # ATTENTION: only works with C_len == 3, C_kernel.shape == (3, 3).
        cooping_competitors = cooping_cells[1, 1:3] if cooping_cells.shape else cooping_cells[2:4, 1]

        # how many cooperators around each competitor?
        C_conv = sp.signal.convolve2d(cooping_cells, self.C_kernel, mode='valid')
        
        
        # Public goods effect.
        # G for Goods.
        # Which are the cells that enjoy the effect of public goods?
        G = (C_conv >= self.C_th)
        
        
        # all cells for which the effect of goods is above threshold is True in G.
        # M for Metabolism
        twocellpos_r, twocellpos_c = sp.arange(rl, rh + 1) % self.N, sp.arange(cl, ch + 1) % self.N
        R_cost_board = self.R_cost * self.B[0, twocellpos_r, twocellpos_c]
        S_cost_board = self.S_cost * self.B[1, twocellpos_r, twocellpos_c]
        C_cost_board = self.C_cost * cooping_competitors
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
        if sp.rand() < self.mutation_rate_r:
            self.B[0, pos[0], pos[1]] = self.B[0, pos[0], pos[1]] ^ True
        if sp.rand() < self.mutation_rate_s:
            self.B[1, pos[0], pos[1]] = self.B[1, pos[0], pos[1]] ^ True
        if sp.rand() < self.mutation_rate_c:
            self.B[2, pos[0], pos[1]] = self.B[2, pos[0], pos[1]] ^ True

    def diffuse(self):
        diffuse_code = r'''
bool tmp_values[3][2][2];

long int genotype_i , row_i , col_i ;
genotype_i = row_i = col_i = 0;

for (row_i = 0; row_i < 2; row_i++)
{
  for (col_i = 0; col_i < 2; col_i++)
  {
    for (genotype_i = 0; genotype_i < 3; genotype_i++)
    {
      tmp_values[genotype_i][row_i % board_size][col_i % board_size] = b[genotype_i, row_i % board_size, col_i % board_size];
    }
  }
}

long int row1 = ((int) row + 1) % board_size;
long int col1 = ((int) col + 1) % board_size;
/*
for (int i; i < 2; i++)
{
  for (int j; i < 2; j++)
  {
    for (int g; g < 3; g++)
    {
      printf("%d", b[g, i, j]);
    }
    printf(" ");
  }
  printf("\n");
}
*/
if (d < 0.5)
{
  for (genotype_i = 0; genotype_i < 3; genotype_i++)
  {
    b[genotype_i, row , col  ] = tmp_values[genotype_i, row , col1]; // anticlockwise index map
    b[genotype_i, row , col1 ] = tmp_values[genotype_i, row1, col1]; // 00 01  01 11
    b[genotype_i, row1, col  ] = tmp_values[genotype_i, row , col ]; // 10 11  00 10
    b[genotype_i, row1, col1 ] = tmp_values[genotype_i, row1, col ];
  }
}
else
{
  for (genotype_i = 0; genotype_i < 3; genotype_i++)
  {
    b[genotype_i, row , col  ] = tmp_values[genotype_i, row1, col ]; // clockwise index map
    b[genotype_i, row , col1 ] = tmp_values[genotype_i, row , col ]; // 00 01  10 00
    b[genotype_i, row1, col  ] = tmp_values[genotype_i, row1, col1]; // 10 11  11 01
    b[genotype_i, row1, col1 ] = tmp_values[genotype_i, row , col1];
  }
}
/*
for (int i; i < 2; i++)
{
  for (int j; i < 2; j++)
  {
    for (int g; g < 3; g++)
    {
      printf("%d", b[g, i, j]);
    }
    printf(" ");
  }
  printf("\n");
}
*/
'''
#        print "yo", self.D, self.N, (self.N ** 2), (self.N ** 2) * self.D, sp.uint64((self.N ** 2) * self.D)
        for i in range(sp.uint64((self.N ** 2) * self.D / 4)):
            direction = sp.random.rand()
            row, col = sp.random.randint(self.N, size=2)
#            print "diffusing:", i
            sp.weave.inline(diffuse_code, ['board_size', 'b', 'row', 'col', 'd'], {'board_size': self.N, 'b': self.B, 'row': row, 'col': col, 'd': direction})
        
    def save_h5(self):
        with h5py.File(self.fname) as ff:
            for key in self.parameters:
                try:
                    print key, type(getattr(self, key))
                    ff[key] = getattr(self, key)
                except:
                    ff[key][...] = getattr(self, key)
    
    def load_h5(self):
        with h5py.File(self.fname) as ff:
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
                nh_genotype_count_of_genotype = sp.sum(nh_genotype_count * genotype_board, dtype=sp.uint64)
                self.samples_nhood[self.sample_count, genotype, nh_genotype] = nh_genotype_count_of_genotype
            self.samples_frequency[self.sample_count, genotype] = genotype_frequency
        self.sample_count += 1

    def nextstep(self):
        print 'generation:', self.step_count
        for i in range(self.N ** 2):
            if i % 100 == 0:
                print 'competition:', i
            self.competition()
        self.diffuse()
        if not self.step_count % self.steps_per_sample:
            self.sample()
        if not self.step_count % self.steps_per_board_sample:
            board_sample_num = self.step_count / self.steps_per_board_sample
            self.samples_board[board_sample_num] = self.B
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

    def display_frequency_timeseries(self):
        for i in range(8):
            pl.plot(sp.arange(self.samples_num), self.samples_frequency[:,i], label=str(i), fillstyle='bottom')
    


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
    signal.signal(signal.SIGINT, handler_maker(a))
    t = time.time()
    every = 30*60
    print "t: %(t)f, steps thus far: %(steps)d" % {'t': t, 'steps': a.step_count}
    steps_a = a.step_count
    while a.step_count <= a.generations:
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            t = time.time()
            a.save_h5()
            steps_delta = a.step_count - steps_a
            steps_a = a.step_count
            eta = 1.0 * delta_t / (steps_delta+1) * (a.steps_final - a.step_count)
            print "t: %(t)f, approx. time to fin: %(eta)f" % {'t': t, 'eta': eta}
            print "steps taken = %(step_count)s, steps since last save = %(steps_delta)s" % {'step_count': a.step_count, 'steps_delta': steps_delta}
            sys.exit(1)

# TODO: Handler of signals.
def handler_maker(a_game):
    def handler(signum, frame):
        print 'Signal handler called with signal', signum
        a_game.save_h5()
        print 'game saved'
        raise
    return handler

if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except IndexError as e:
        print "Usage: %(me)s [datafilename]" % { 'me' : sys.argv[0] }
        raise
    a = game_of_strife(fname = fname, N=10, generations=1)
    if os.path.exists(fname):
        a.load_h5()
        go(a)
    else:
        a = game_of_strife(fname = fname, N=100, generations=10000)
        a.save_h5()
        go(a)
    a.save_h5()
    sys.exit(0)
