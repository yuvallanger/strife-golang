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
import ConfigParser

labels = ['Ignorant (csr)', 'Voyeur (csR)', 'Liar (cSr)', 'Lame (cSR)',
          'Blunt (Csr)', 'Shy (CsR)', 'Vain (CSr)', 'Honest (CSR)']

class gameOfStrife:
    def __init__(self, config=None):
        if config == None:
            config = default_config

        d = {}

        ########
        # filname
        d['fname'] = config['data_filename']

        #########
        # model parameters
        #########
        
        #######
        # Cost of gene expression
        d['S_cost'] = sp.int64(config['S_cost'])
        d['R_cost'] = sp.int64(config['R_cost'])
        d['C_cost'] = sp.int64(config['C_cost'])
        d['metabolic_baseline'] = sp.int64(config['metabolic_baseline'])  # B for Baseline, basal, "basic metabolic burden"

        #####
        # benefit from cooperation. "reward factor" in the article.
        d['benefit'] = sp.float64(config['benefit'])

        ######
        # mutation per generation
        d['mutation_rate_r'] = sp.float64(config['mutation_rate_r'])
        d['mutation_rate_s'] = sp.float64(config['mutation_rate_s'])
        d['mutation_rate_c'] = sp.float64(config['mutation_rate_c'])

        ## neighbours effects' thresholds
        d['S_th'] = sp.int64(config['S_th'])
        # quorum threshold
        d['C_th'] = sp.int64(config['C_th'])
        # Cooperation threshold. Above it, public goods makes a difference.
        
        ## Probability of each single diffusion operation
        d['diffusion_amount'] = sp.float64(config['diffusion_amount'])
        
        #######
        # radius of Signal or Cooperation effects.
        d['S_rad'] = sp.int64(config['S_rad'])
        d['C_rad'] = sp.int64(config['C_rad'])
       
        d['generations'] = sp.int64(config['generations'])

        d['N'] = sp.int64(config['board_size'])
      
        #####
        # settings
        #####
        
        d['NEIGHBOUR_REL_POS'] = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        #######
        # we'll increase step_count by one every time two cells compete.
        d['step_count'] = sp.int64(0)


        d['steps_final'] = sp.int64(d['generations'] * d['N']**2)

        # diameter of the convolution matrix
        diameter = lambda x: sp.int64(2 * x + 1)
        S_len = diameter(d['S_rad'])
        C_len = diameter(d['C_rad'])

        # the convolution matrix used to count neighbours
        d['S_kernel'] = sp.ones((S_len, S_len))
        d['C_kernel'] = sp.ones((C_len, C_len))

        # A cell can be Signalling and/or Receptive and/or Cooperative
        R = sp.rand(d['N'], d['N']) > config['initial_receptives_amount']
        S = sp.rand(d['N'], d['N']) > config['initial_signallers_amount']
        C = sp.rand(d['N'], d['N']) > config['initial_cooperators_amount']
        d['B'] = sp.array([R, S, C]).transpose((1,2,0))
        assert d['N'] == d['B'].shape[0] and d['N'] == d['B'].shape[1], 'B.shape: {0}\nN: {1}\nWanted: {2}'.format(d['B'].shape, d['N'], (d['N'], d['N'], 3))
        
        d['genotype_num'] = sp.int64(8)
        
        ## data sampling
        # we will take a frequency sample some number of times per generation
        d['steps_per_gen'] = sp.int64(d['N'] ** 2)
        d['samples_per_gen'] = sp.int64(1)
        d['samples_num'] = sp.int64(d['samples_per_gen'] * d['generations'])
        d['samples_board_num'] = sp.int64(d['samples_num'] // 10)
        d['steps_per_sample'] = sp.int64(sp.floor(1.0 * d['steps_per_gen'] // d['samples_per_gen']), dtype=sp.int64)
        d['steps_per_board_sample'] = sp.int64(10 * d['steps_per_sample'])
        d['sample_count'] = sp.int64(0)
        # We want to know the frequency of each genotype per generation
        d['samples_frequency'] = sp.empty((d['samples_num'], d['genotype_num']), dtype='int32')
        d['samples_nhood'] = sp.empty((d['samples_num'], d['genotype_num'], d['genotype_num']), dtype=sp.int64)
        d['samples_board'] = sp.empty((d['samples_board_num'], d['N'], d['N'], 3), dtype=sp.int64)
       
        if not config == None:
            self.__init_unpack_parameters__(d)
    
    def __init_unpack_parameters__(self, d):

        """
        __init_unpack_parameters__(self, d)

        Introdueces the key-value pairs of dictionary "d" as attributes of self.
        """

        self.parameters = set()
        for key, val in d.items():
            setattr(self, key, val)
            self.parameters.add(key)

    ######
    ## functions
    ######    
    
#    @profile
    def competition(self, c_pos_1, c_pos_2, p_pair):

        """
        competition(self, cell_pos_1, cell_pos_2, p_pair) -> (int, int)

        Decides which of the two positions wins.

        Coordinates are a numpy array of shape=(2,) and an integer dtype.
        Takes two adjacent position coordinates on the board and each one's TODO: what's the name of such probability values?
        p_pair: an array of two uniform distribution over [0, 1).
        """

        p1, p2 = p_pair
        assert p1.shape=(2,) & p2.shape=(2,), 'p1 ({0}) and p2 ({1}) need to be of shape (2,)'.format(p1, p2)
        assert p1.dtype.kind == 'i' & p2.dtype.kind == 'i', 'p1 ({0}) and p2 ({1}) need to be of integer dtype'.format(p1, p2)
        assert ((0 <= p1) & (p1 < 1) & (0 <= p2) & (p2 < 1)).all(), 'p1 ({0}) and p2 ({1}) need to be over [0, 1)'.format(p1, p2)

        # c_pos_2's coordinates in a torus:
        c_pos_2t = c_pos_2 % self.board_edge_length
        assert (0 <= c_pos_1[0]) and \
               (0 <= c_pos_1[1]) and \
               (0 <= c_pos_2t[0]) and \
               (0 <= c_pos_2t[1]) and \
               (c_pos_1[0] < self.board_edge_length) and \
               (c_pos_1[1] < self.board_edge_length) and \
               (c_pos_2t[0] < self.board_edge_length) and \
               (c_pos_2t[1] < self.board_edge_length), 'c_pos_1: {0}\nc_pos_2t: {1}'.format(c_pos_1, c_pos_2t)

        # two identical cells competing will result in two identical cells,
        # so we will return now with no further calculation of this competition.
        if (self.board[c_pos_1[0], c_pos_1[1]] == self.board[c_pos_2t[0], c_pos_2t[1]]).all():
            return (c_pos_2t, c_pos_1)

        ## We will optimize by taking a sub array from each genotype array around the competitors.

        # rl, ch - row low, col high
#        twosort = lambda x, y: (x, y) if x < y else (y, x)
        rl, rh = (c_pos_1[0], c_pos_2[0]) if c_pos_1[0] < c_pos_2[0] else (c_pos_2[0], c_pos_1[0])
        cl, ch = (c_pos_1[1], c_pos_2[1]) if c_pos_1[1] < c_pos_2[1] else (c_pos_2[1], c_pos_1[1])

        # I'll use this to keep the arrays as 2d (ndim=2)
        assert rl == rh or cl == ch, 'rl: {0}\nrh: {1}\ncl: {2}\n ch: {3}'.format(rl, rh, cl, ch)
        
        if rl == rh:
            shape = (1,2)
        elif cl == ch:
            shape = (2,1)
        
        # For signallers, we take both S_rad and C_rad around our competitors because
        # signallers affect CR[Ss] cells which, with their public goods, affect our competitors
        s_row_range = sp.arange(rl - self.S_rad - self.C_rad, rh + self.S_rad + self.C_rad + 1) % self.board_edge_length
        s_col_range = sp.arange(cl - self.S_rad - self.C_rad, ch + self.S_rad + self.C_rad + 1) % self.board_edge_length
        rc_row_range = sp.arange(rl - self.C_rad, rh + self.C_rad + 1) % self.board_edge_length
        rc_col_range = sp.arange(cl - self.C_rad, ch + self.C_rad + 1) % self.board_edge_length

        assert self.S_rad.dtype.kind == sp.int8(1).dtype.kind, 'Got {0}, wanted {1}'.format(self.S_rad.dtype.kind, sp.int_(1).dtype.kind)
        assert self.C_rad.dtype.kind == sp.int8(1).dtype.kind, 'Got {0}, wanted {1}'.format(self.C_rad.dtype.kind, sp.int_(1).dtype.kind)
        assert s_row_range.shape == (shape[0]-1 + 2 * self.S_rad + 2 * self.C_rad + 1, ), '''s_row_range: {0}
rl: {1}; rh: {2}
wanted: {3}'''.format(s_row_range, rl, rh, (shape[0]-1 + 2 * self.S_rad + 2 * self.C_rad + 1, ))
        assert s_col_range.shape == (shape[1]-1 + 2 * self.S_rad + 2 * self.C_rad + 1, ), '''s_col_range: {0}
cl: {1}; ch: {2}
wanted: {3}'''.format(s_col_range, cl, ch, (shape[1]-1 + 2 * self.S_rad + 2 * self.C_rad + 1, ))
        assert rc_row_range.shape == (shape[0]-1  + 2 * self.C_rad + 1, ), '''rc_row_range: {0}
rl: {1}; rh: {2}
wanted: {3}'''.format(rc_row_range, rl, rh, (shape[0]-1 + 2 * self.C_rad + 1, ))
        assert rc_col_range.shape == (shape[1]-1 + 2 * self.C_rad + 1, ), '''rc_col_range: {0}
cl: {1}; ch: {2}
wanted: {3}'''.format(rc_col_range, cl, ch, (shape[1]-1 + 2 * self.C_rad + 1, ))

        R_sub = self.board[rc_row_range, :, 0][:, rc_col_range]
        S_sub = self.board[s_row_range, :, 1][:, s_col_range]
        C_sub = self.board[rc_row_range, :, 2][:, rc_col_range]

        assert R_sub.shape == tuple(sp.array(shape)+2), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+2)
        assert S_sub.shape == tuple(sp.array(shape)+4), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+4)
        assert C_sub.shape == tuple(sp.array(shape)+2), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+2)

        # we count how many signallers are within each cell's neighbourhood
        S_conv = sp.signal.convolve2d(S_sub, self.S_kernel, mode='valid')

        assert_ndim(S_conv, 2)

        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
        cooping_cells = ((C_sub & R_sub) & (S_conv >= self.S_th)) | (C_sub & (R_sub ^ True))

        assert_ndim(cooping_cells, 2) # TODO: Continue putting loads of assert_ndim()s
        assert (cooping_cells.shape == (3,4) and shape == (1,2)) or (cooping_cells.shape == (4,3) and shape == (2,1)), '''cooping_cells.shape: {0}
shape: {1}'''.format(cooping_cells.shape, shape)

        # ATTENTION: only works with C_len == 3, C_kernel.shape == (3, 3).
        if cooping_cells.shape == (4, 3):
            cooping_competitors = cooping_cells[1:3, 1].reshape(shape)
        elif cooping_cells.shape == (3, 4):
            cooping_competitors = cooping_cells[1, 1:3].reshape(shape)
        else:
            raise 'blah'

        assert cooping_competitors.shape == shape, 'cooping_competitors: {0}\nshape: {1}\nWanted shape: (1,2) or (2,1)'.format(cooping_competitors.shape, shape)

        # how many cooperators around each competitor?
        C_conv = sp.signal.convolve2d(cooping_cells, self.C_kernel, mode='valid')
        assert C_conv.shape == shape, '''C_conv.shape: {0}
shape: {1}'''.format(C_conv.shape, shape)
        
        # Public goods effect.
        # G for Goods.
        # Which are the cells that enjoy the effect of public goods?
        G = (C_conv >= self.C_th)

        assert G.shape == shape, 'G.shape: {0}\nshape: {1}'.format(G.shape, shape)
        
        # all cells for which the effect of goods is above threshold is True in G.
        # M for Metabolism
        twocellpos_r = sp.arange(rl, rh + 1) % self.board_edge_length
        twocellpos_c = sp.arange(cl, ch + 1) % self.board_edge_length

        assert (twocellpos_r.shape == (2,) and twocellpos_c.shape == (1,)) or (twocellpos_r.shape == (1,) and twocellpos_c.shape == (2,)), 'twocellpos_r.shape: {0}\ntwocellpos_c.shape: {1}\nshape: {2}'.format(twocellpos_r.shape, twocellpos_c.shape, shape)

        R_cost_board = self.R_cost * self.board[twocellpos_r, twocellpos_c, 0].reshape(shape)
        S_cost_board = self.S_cost * self.board[twocellpos_r, twocellpos_c, 1].reshape(shape)
        C_cost_board = self.C_cost * cooping_competitors

        assert R_cost_board.shape == shape, 'R_cost_board: {0}\nWanted ndim: 1'.format(R_cost_board)
        assert S_cost_board.shape == shape, 'S_cost_board: {0}\nWanted ndim: 1'.format(S_cost_board)
        assert C_cost_board.shape == shape, 'C_cost_board: {0}\nWanted ndim: 1'.format(C_cost_board)

        Total_cost_board = S_cost_board + R_cost_board + C_cost_board + self.metabolic_baseline

        assert_ndim(Total_cost_board, 2)

        M = G * (1 - self.benefit) * Total_cost_board
        assert_ndim(M, 2)
        # all false in G don't benefit from public goods (G^True flips values)
        M += (G^True) *  Total_cost_board
        assert_ndim(M, 2)
        M = self.metabolic_baseline / M
        assert_ndim(M, 2)
        score1 = p1 * M.item(0) # score1 is the first position's score
        score2 = p2 * M.item(1) # score2 is the second position's score
        if shape == (2, 1):
            assert M.shape == (2, 1), 'M.shape == {0}\nWanted == {1}'.format(M.shape, (2, 1))
            if c_pos_2[0] > c_pos_1[0]:
                # their position is like this:
                # 2
                # 1
                if score1 > score2:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
                else:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1
                # 2
                if score1 > score2:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
                else:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
        else:
            assert M.shape == (1, 2), 'M.shape == {0}\nWanted == {1}'.format(M.shape, (1, 2))
            if c_pos_2[1] < c_pos_1[1]:
                # their position is like this:
                # 2 1
                if score1 > score2:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
                else:
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1 2
                if score1 > score2:
                    #print "comp 1 wins"
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
                else:
                    #print "comp 2 wins"
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
    
    def copycell(self, orig, dest):

        """
        copycell(self, orig, dest) -> NoneType

        Copies the contents of self.board at coordinates of position "orig" into the position of coordinates "dest".
        Coordinates are a numpy array of shape=(2,) and an integer dtype.
        """

        assert orig.shape == (2,) and dest.shape == (2,), 'orig.shape: {0}\ndest.shape: {1}'.format(orig.shape, dest.shape)
        self.board[dest[0], dest[1]] = self.board[orig[0], orig[1]]

    def mutate(self, pos):

        """
        mutate(self, pos) -> NoneType

        For each value of self.board at position "pos", change its value at probability of self.mutation_rate_[r/s/c].
        """

        if sp.rand() < self.mutation_rate_r:
            self.board[pos[0], pos[1], 0] = self.board[pos[0], pos[1], 0] ^ True
        if sp.rand() < self.mutation_rate_s:
            self.board[pos[0], pos[1], 1] = self.board[pos[0], pos[1], 1] ^ True
        if sp.rand() < self.mutation_rate_c:
            self.board[pos[0], pos[1], 2] = self.board[pos[0], pos[1], 2] ^ True

    def diffuse(self, direction, position):

        """
        diffuse(self, direction, position) -> NoneType

        Turns a 2 by 2 sub-array of self.board by 90 degrees.
        Direction:
            True - turn 
        with its lowest valued coordinate (upper left) in the "position" coordinate value.
        """
        
        diffuse_code = r'''
//bool temp_values[3][2][2];

long int genotype_i , row_i , col_i ;
genotype_i = row_i = col_i = 0;

for (row_i = 0; row_i < 2; row_i++)
{
  for (col_i = 0; col_i < 2; col_i++)
  {
    for (genotype_i = 0; genotype_i < 3; genotype_i++)
    {
      temp_values[row_i % board_size][col_i % board_size][genotype_i] = b[row_i % board_size, col_i % board_size, genotype_i];
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
      printf("%d", b[i, j, g]);
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
    b[ row , col , genotype_i ] = temp_values[ row , col1, genotype_i ]; // anticlockwise index map
    b[ row , col1, genotype_i ] = temp_values[ row1, col1, genotype_i ]; // 00 01 -> 01 11
    b[ row1, col , genotype_i ] = temp_values[ row , col , genotype_i ]; // 10 11 -> 00 10
    b[ row1, col1, genotype_i ] = temp_values[ row1, col , genotype_i ];
  }
}
else
{
  for (genotype_i = 0; genotype_i < 3; genotype_i++)
  {
    b[ row , col , genotype_i ] = temp_values[ genotype_i, row1, col  ]; // clockwise index map
    b[ row , col1, genotype_i ] = temp_values[ genotype_i, row , col  ]; // 00 01 -> 10 00
    b[ row1, col , genotype_i ] = temp_values[ genotype_i, row1, col1 ]; // 10 11 -> 11 01
    b[ row1, col1, genotype_i ] = temp_values[ genotype_i, row , col1 ];
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
        row, col = position
        temp_values = sp.empty((3,2,2), dtype=sp.bool_)
        
        sp.weave.inline(diffuse_code, ['board_size',
                                       'board',
                                       'row',
                                       'col',
                                       'direction',
                                       'temp_valuee'],
                                       {'board_size': self.board_edge_length,
                                        'board': self.board,
                                        'row': row,
                                        'col': col,
                                        'direction': direction
                                        'temp_value': temp_value })
        
    def save_h5(self):

        """
        Saves the attributes of self whose names show up as keys in self.parameters.
        """

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
        joint_board = self.board[:, :, 0] + 2 * self.board[:, :, 1] + 4 * self.board[:, :, 2]
        for genotype in range(8):
            genotype_board = joint_board == genotype
            genotype_frequency = sp.sum(genotype_board)
            # neighbours_genotype
            for nh_genotype in range(8):
                nh_board = joint_board == nh_genotype
                nh_genotype_count = sp.signal.convolve2d(nh_board, sp.ones((3,3)), mode='same', boundary='wrap')
                nh_genotype_count_of_genotype = sp.sum(nh_genotype_count * genotype_board, dtype=sp.int64)
                self.samples_nhood[self.sample_count, genotype, nh_genotype] = nh_genotype_count_of_genotype
            self.samples_frequency[self.sample_count, genotype] = genotype_frequency
        self.sample_count += 1

    def nextstep(self):
        print 'generation:', self.step_count
        for i in range(self.board_edge_length ** 2):
            print 'competition:', i
            ##
            # Draw two adjacent positions.
            # We'll use relative positions to compute exact positions of 2nd competitor cell
            pos1 = sp.random.randint(self.board_edge_length, size=2)
            pos2 = pos1 + self.NEIGHBOUR_REL_POS[sp.random.randint(4)]
            p_pair = sp.rand(2)
            winner_pos, loser_pos = self.competition(pos1, pos2, p_pair)
            self.copycell(winner_pos, loser_pos)
            self.mutate(loser_pos)
        for i in range(sp.int64((self.board_edge_length ** 2) * self.diffusion_amount // 4)):
            print 'diffusion: {0}'.format(i)
            direction = sp.rand()
            position = sp.random.randint(self.board_edge_length, size=2)
            self.diffuse(direction, position)
        if not self.step_count % self.steps_per_sample:
            self.sample()
        if not self.step_count % self.steps_per_board_sample:
            board_sample_num = self.step_count // self.steps_per_board_sample
            self.samples_board[board_sample_num] = self.board
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
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    print "t: {0:f}, steps thus far: {1:d}".format(t, a.step_count)
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
            print "t: {0:f}, approx. time to fin: {1:f}".format(t, eta)
            print "steps taken = {0}, steps since last save = {1}".format(a.step_count, steps_delta)
            sys.exit(1)

# TODO: Handler of signals.
def handler_maker(a_game):
    def handler(signum, frame):
        print 'Signal handler called with signal', signum
        a_game.save_h5()
        print 'game saved'
        raise
    return handler

def assert_ndim(arr, nd):
    assert arr.ndim == nd, 'Wrong number of dimensions.\nExists {0} but {1} is wanted.'.format(arr.ndim, nd)

def assert_shape(arr, shape):
    assert arr.shape == shape, 'Wrong shape\nExists {0} but {1} is wanted.'.format(arr.shape, shape)

default_config = {
        'S_cost':             1, # Metabolic cost of signalling
        'R_Metabolic cost':   3, # Metabolic cost of having a receptor
        'C_Metabolic cost':  30, # Metabolic cost of being cooperative
        'metabolic_baseline':           100, # Basal metabolic cost
        'benefit':          0.9, # The fraction reduced, out of total metabolic cost, when
                                 #    public goods reach threshold of benefit.
        # Likelihoods of switch (on to off and vica versa) for each gene per cell per generation.
        'mutation_rate_r': 1e-4,
        'mutation_rate_s': 1e-4,
        'mutation_rate_c': 1e-4,
        'S_th':               3, # Amount of signal needed for a receptive and cooperative cell to
                                 #    start producing public goods.
        'C_th':               3, # Amount of public goods needed for metabolic benefit.
        'diffusion_amount': 0.5, # Average fraction of cells out of total cells on the board (board_size**2)
                                 #    which will be involved
        'board_size':        10, # The length of the board. The board will be of shape (board_size, board_size).
        'generations':       10, # Each generation involves, on average, all cells in a competition, meaning
                                 #    board_size**2 competitions.
        'S_rad':              1, # Radius of the signal's effect.
        'C_rad':              1,
        'samples_per_gen':    1,
        'initial_receptives_amount': 0,
        'initial_signallers_amount': 0,
        'initial_cooperators_amount': 0,
        'data_filename': 'strife.h5' }

def load_config(config_filename):
    """
    Takes a string holding filename and returns a dictionary with all the configuration values.
    """
    our_config = default_config

    config = ConfigParse.SafeConfigParser()
    config.read(config_filename)

    for key, val in our_config:
        our_config[key] = config.get('Config', key)

    return our_config

