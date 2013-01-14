# -*- coding: utf-8 -*-
# cython: profile=True

"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import os
import time
import scipy as sp
#cimport scipy as sp
import scipy.signal
import scipy.weave
import numpy as np
cimport numpy as np
#import pygame
import pylab as pl
#import timeit
import sys
import ConfigParser
cimport cython

cdef int SIGNAL, RECEPTOR, COOPERATION
RECEPTOR = 0
SIGNAL = 1
COOPERATION = 2

cdef double S_cost
cdef double R_cost
cdef double C_cost
cdef double baseline_cost
cdef double benefit
cdef double mutation_rate_r
cdef double mutation_rate_s
cdef double mutation_rate_c
cdef int S_th
cdef int G_th
cdef double diffusion_amount
cdef int S_rad
cdef int G_rad
cdef int generations
cdef int board_size
cdef int step_count
cdef int steps_final
cdef np.ndarray board
cdef char[:,:,:] board_memview
cdef int genotype_num
cdef int steps_per_gen
cdef int samples_per_gen
cdef int samples_num
cdef int samples_board_num
cdef int steps_per_sample
cdef int steps_per_board_sample
cdef int sample_count
cdef np.ndarray samples_frequency
cdef np.ndarray samples_nhood
cdef np.ndarray samples_board

def init( config):
    global S_cost, R_cost, C_cost, baseline_cost, benefit, mutation_rate_r, mutation_rate_s, mutation_rate_c, S_th, G_th, diffusion_amount, S_rad, G_rad, generations, board_size, step_count, steps_final, board, board_memview, genotype_num, steps_per_gen, samples_per_gen, samples_num, samples_board_num, steps_per_sample, steps_per_board_sample, sample_count, samples_frequency, samples_nhood, samples_board
    #########
    # model parameters,
    #########
    
    #######
    # Cost of gene expression
    S_cost = np.int32(config['S_cost'])
    R_cost = np.int32(config['R_cost'])
    C_cost = np.int32(config['C_cost'])
    baseline_cost = np.int32(config['baseline_cost'])  # B for Baseline, basal, "basic metabolic burden"

    #####
    # benefit from cooperation. "reward factor" in the article.
    benefit = np.float64(config['benefit'])

    ######
    # mutation per generation
    mutation_rate_r = np.float64(config['mutation_rate_r'])
    mutation_rate_s = np.float64(config['mutation_rate_s'])
    mutation_rate_c = np.float64(config['mutation_rate_c'])

    ## neighbours effects' thresholds
    S_th = np.int32(config['S_th'])
    # quorum threshold
    G_th = np.int32(config['G_th'])
    # Cooperation threshold. Above it, public goods makes a difference.
    
    ## Probability of each single diffusion operation
    diffusion_amount = np.float64(config['diffusion_amount'])
    
    #######
    # radius of Signal or Cooperation effects.
    S_rad = np.int32(config['S_rad'])
    G_rad = np.int32(config['G_rad'])
   
    generations = np.int32(config['generations'])

    board_size = np.int32(config['board_size'])
  
    #####
    # settings
    #####
    
    ## time keeping
    # number of generations the simulation will run
    # each generation is defined as the average number of steps for which
    # each cell on the board was in a competition once since last
    # generation (?)

    #######
    # we'll increase step_count by one every time two cells compete.
    step_count = np.int32(0)


    steps_final = np.int32(generations * board_size**2)

    # A cell can be Signalling and/or Receptive and/or Cooperative
    R = np.random.rand(board_size, board_size) < np.float64(config['initial_receptives_amount'])
    S = np.random.rand(board_size, board_size) < np.float64(config['initial_signallers_amount'])
    C = np.random.rand(board_size, board_size) < np.float64(config['initial_cooperators_amount'])
    R = np.int8(R)
    S = np.int8(S)
    C = np.int8(C)
    board = np.array([R, S, C]).transpose((1,2,0))
    board_memview = board

    # TODO do whole board sum
    # S_sum = np.empty((board_size, board_size), 
    
    genotype_num = np.int32(8)
    
    ## data sampling
    # we will take a frequency sample some number of times per generation
    steps_per_gen = np.int32(board_size ** 2)
    samples_per_gen = np.int32(1)
    samples_num = np.int32(samples_per_gen * generations)
    samples_board_num = np.int32(samples_num // 10)
    steps_per_sample = np.int32(np.floor(1.0 * steps_per_gen // samples_per_gen))
    steps_per_board_sample = np.int32(10 * steps_per_sample)
    sample_count = np.int32(0)
    # We want to know the frequency of each genotype per generation
    samples_frequency = np.empty((samples_num, genotype_num), dtype=np.int32)
    samples_nhood = np.empty((samples_num, genotype_num, genotype_num), dtype=np.int32)
    samples_board = np.empty((samples_board_num, board_size, board_size, 3), dtype=np.int8)

######
## functions
######
cdef int signal_count(int row, int col, int signal):
    cdef int signal_sum, s_row_i, s_col_i
    signal_sum = 0
    for s_row_i in range(row - self.S_rad, row + self.S_rad + 1):
        for s_col_i in range(col - self.S_rad, col + self.S_rad + 1):
            if board_memview[s_row_i, s_col_i, SIGNAL] == signal:
                signal_sum += 1
    return signal_sum

cdef int goods_count(int pos_row, int pos_col):
    '''
    Counts the number of public goods at a certain position
    '''
    cdef int s_row_i, s_col_i, g_row_i, g_col_i, goods_sum
    
#     has fitness effect?:
#         surrounded by enough goods producing cells
#     produces goods?:
#         no receptor and cooperator
#         or
#         surrounded by enough signal producing cells and
#             is receptive and
#             is cooperator

    for g_row_i in range(pos_row - G_rad,
                         pos_row + G_rad + 1):
        for g_col_i in range(pos_col - G_rad,
                             pos_col + G_rad + 1):
            # Able to produce public goods
            if board_memview[g_row_i % self.board_size,
                     g_col_i % self.board_size,
                     COOPERATION] == 1:
                # Isn't receptive
                if board_memview[g_row_i % self.board_size,
                         g_col_i % self.board_size,
                         RECEPTOR] == 0:
                    goods_sum += 1
                # Receptive and signal reaches signal threshold
                elif self.signal_count(g_row_i, g_col_i, 1) >= self.S_th:
                    goods_sum += 1
    return goods_sum

cdef double fitness(self, pos_row, pos_col): 
    cdef double result
    cdef int pos_row_t, pos_col_t
    pos_row_t, pos_col_t = (pos_row % self.board_size,
                            pos_row % self.board_size)

    result = self.benefit * (self.goods_count(pos_row_t, pos_col_t) >= self.G_th)
    result = 1 - result
    result = (self.baseline_cost + result) / result
    result = result

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cdef competition(self,
                  int c_pos_1_row,
                  int c_pos_1_col,
                  int c_pos_2_row,
                  int c_pos_2_col,
                  double p_1,
                  double p_2):
    """
    competition(self, cell_pos_1, cell_pos_2, p_pair) -> (int, int)

    Decides which of the two positions wins.

    Coordinates are a numpy array of shape=(2,) and an integer dtype.
    Takes two adjacent position coordinates on the board and each one's TODO: what's the name of such probability values?
    p_pair: an array of two uniform distribution over [0, 1).
    """
#    cdef double score_1, score_2
#    cdef int c_pos_2t_row, c_pos_2t_col
#    c_pos_2t_row = c_pos_2_row % self.board_size 
#    c_pos_2t_col = c_pos_2_col % self.board_size
#
#    score_1 = p_1 * self.fitness(c_pos_1_row,
#                                 c_pos_1_col)
#    score_2 = p_2 * self.fitness(c_pos_2t_row,
#                                 c_pos_2t_col)
#
#    return 1 if score_1 > score_2 else 2
    return 1
    
cdef copycell(self,
             int orig_row, int orig_col,
             int dest_row, int dest_col):
    """
    copycell(self, orig_row, orig_col, dest_row, dest_col) -> NoneType

    Copies the contents of self.board at coordinates of position "orig" into the position of coordinates "dest".
    Coordinates are a numpy array of shape=(2,) and an integer dtype.
    """
    cdef int i, orig_row_t, orig_col_t, dest_row_t, dest_col_t
    orig_row_t = orig_row % self.board_size
    orig_col_t = orig_col % self.board_size
    dest_row_t = dest_row % self.board_size
    dest_col_t = dest_col % self.board_size

    for i in range(3):
        board_memview[dest_row_t, dest_col_t, i] = board_memview[orig_row_t, orig_col_t, i]

def mutate(self, pos_row, pos_col, p_r, p_s, p_c):
    """
    mutate(self, pos) -> NoneType

    For each value of self.board at position "pos", change its value at probability of self.mutation_rate_[r/s/c].
    """
    cdef pos_row_t, pos_col_t
    pos_row_t = pos_row % self.board_size
    pos_col_t = pos_col % self.board_size
    if p_r < self.mutation_rate_r:
        board_memview[pos_row_t, pos_col_t, 0] = 0 if board_memview[pos_row_t, pos_col_t, 0] else 1
    if p_s < self.mutation_rate_s:
        board_memview[pos_row_t, pos_col_t, 1] = 0 if board_memview[pos_row_t, pos_col_t, 1] else 1
    if p_c < self.mutation_rate_c:
        board_memview[pos_row_t, pos_col_t, 2] = 0 if board_memview[pos_row_t, pos_col_t, 2] else 1

def get_state(self):
    """
    The state of the simulation in a key, val dictionary.
    """
    return {'S_cost': S_cost,
            'R_cost': R_cost,
            'C_cost': C_cost,
            'baseline_cost': baseline_cost,
            'benefit': benefit,
            'mutation_rate_r': mutation_rate_r,
            'mutation_rate_s': mutation_rate_s,
            'mutation_rate_c': mutation_rate_c,
            'S_th': S_th,
            'G_th': G_th,
            'diffusion_amount': diffusion_amount,
            'S_rad': S_rad,
            'G_rad': G_rad,
            'generations': generations,
            'board_size': board_size,
            'step_count': step_count,
            'steps_final': steps_final,
            'board': board,
            'genotype_num': genotype_num,
            'steps_per_gen': steps_per_gen,
            'samples_per_gen': samples_per_gen,
            'samples_num': samples_num,
            'samples_board_num': samples_board_num,
            'steps_per_sample': steps_per_sample,
            'steps_per_board_sample': steps_per_board_sample,
            'sample_count': sample_count,
            'samples_frequency': samples_frequency,
            'samples_nhood': samples_nhood,
            'samples_board': samples_board}

def set_state( data):
    """
    Sets the state of the simulation with a key, val dictionary holding the data.
    """
    global S_cost, R_cost, C_cost, baseline_cost, benefit, mutation_rate_r, mutation_rate_s, mutation_rate_c, S_th, G_th, diffusion_amount, S_rad, G_rad, generations, board_size, step_count, steps_final, board, board_memview, genotype_num, steps_per_gen, samples_per_gen, samples_num, samples_board_num, steps_per_sample, steps_per_board_sample, sample_count, samples_frequency, samples_nhood, samples_board
    S_cost = data['S_cost']
    R_cost = data['R_cost']
    C_cost = data['C_cost']
    baseline_cost = data['baseline_cost']
    benefit = data['benefit']
    mutation_rate_r = data['mutation_rate_r']
    mutation_rate_s = data['mutation_rate_s']
    mutation_rate_c = data['mutation_rate_c']
    S_th = data['S_th']
    G_th = data['G_th']
    diffusion_amount = data['diffusion_amount']
    S_rad = data['S_rad']
    G_rad = data['G_rad']
    generations = data['generations']
    board_size = data['board_size']
    step_count = data['step_count']
    steps_final = data['steps_final']
    board = data['board']
    genotype_num = data['genotype_num']
    steps_per_gen = data['steps_per_gen']
    samples_per_gen = data['samples_per_gen']
    samples_num = data['samples_num']
    samples_board_num = data['samples_board_num']
    steps_per_sample = data['steps_per_sample']
    steps_per_board_sample = data['steps_per_board_sample']
    sample_count = data['sample_count']
    samples_frequency = data['samples_frequency']
    samples_nhood = data['samples_nhood']
    samples_board = data['samples_board']

def sample():
    global samples_nhood, samples_frequency, sample_count
    joint_board = board[:, :, 0] + 2 * board[:, :, 1] + 4 * board[:, :, 2]
    for gene in range(8):
        gene_board = joint_board == gene
        gene_frequency = np.sum(gene_board)
        # neighbours_genotype
        for nh_gene in range(8):
            nh_board = joint_board == nh_gene
            nh_gene_count = sp.signal.convolve2d(nh_board, np.ones((3,3)), mode='same', boundary='wrap')
            nh_gene_count_of_gene = np.sum(nh_gene_count * gene_board, dtype=np.int32)
            samples_nhood[sample_count, gene, nh_gene] = nh_gene_count_of_gene
        samples_frequency[sample_count, gene] = gene_frequency
    sample_count += 1

cdef int rel_pos_find(self, int i, int axis):
    if i == 0:
        if axis == 0:
            return -1
        else:
            return 0
    elif i == 0:
        if axis == 0:
            return 0
        else:
            return -1
    elif i == 0:
        if axis == 0:
            return 1
        else:
            return 0
    else:
        if axis == 0:
            return 0
        else:
            return 1

cdef same(self, a_r, a_c, b_r, b_c):
    cdef i
    for i in range(3):
        if not (self.board_memview[a_r, a_c, i] == self.board_memview[b_r, b_c, i]):
            return 0
    return 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef nextstep(self):
    cdef int i
    #cdef char[:,:,:] board_memview = self.board
    cdef int rel_pos_i, pos_1_row, pos_1_col, pos_2_row, pos_2_col, winner, pos_2t_row, pos_2t_col
    cdef double p_1, p_2, p_r, p_s, p_c
    cdef np.ndarray[np.int_t, ndim=2] pos_1s
    cdef int[:,:] pos_1s_memview
    cdef np.ndarray[np.int_t, ndim=1] rel_pos_s
    cdef int[:] rel_pos_s_memview
    cdef np.ndarray[np.double_t, ndim=2] p_pairs
    cdef double[:,:] p_pairs_memview
    cdef np.ndarray[np.double_t, ndim=2] p_muts
    cdef double[:,:] p_muts_memview

    pos_1s = np.random.randint(self.board_size, size=(self.board_size ** 2, 2))
    rel_pos_s = np.random.randint(4, size=(self.board_size ** 2))
    p_pairs = np.random.rand(self.board_size ** 2,  2)
    p_muts = np.random.rand(self.board_size ** 2, 3)
    for i from 0 <= i < (self.board_size ** 2):
        pos_2_row = pos_1s[i, 0] + self.rel_pos_find(rel_pos_s[i], 0)
        pos_2_col = pos_1s[i, 1] + self.rel_pos_find(rel_pos_s[i], 1)
        pos_2t_row = pos_2_row % self.board_size
        pos_2t_col = pos_2_col % self.board_size

        if self.same(pos_1s[i,0],
                     pos_1s[i,1],
                     pos_2t_row,
                     pos_2t_col):
            self.mutate(pos_1s[i, 0],
                        pos_1s[i, 1],
                        p_muts[i, 0],
                        p_muts[i, 1],
                        p_muts[i, 2])
            continue

        winner = self.competition(pos_1s[i, 0],
                                  pos_1s[i, 1],
                                  pos_2_row,
                                  pos_2_col,
                                  p_pairs[i, 0],
                                  p_pairs[i, 1])
        if winner == 1:
            self.copycell(pos_1s[i, 0],
                          pos_1s[i, 1],
                          pos_2_row,
                          pos_2_col)
            self.mutate(pos_2_row,
                        pos_2_col,
                        p_muts[i, 0],
                        p_muts[i, 1],
                        p_muts[i, 2])
        else:
            self.copycell(pos_2_row,
                          pos_2_col,
                          pos_1s[i, 0],
                          pos_1s[i, 1])
            self.mutate(pos_1s[i, 0],
                        pos_1s[i, 1],
                        p_muts[i, 0],
                        p_muts[i, 1],
                        p_muts[i, 2])
#   for i in range(np.int32((self.board_size ** 2) * self.diffusion_amount // 4)):
#       print 'diffusion: {0}'.format(i)
#       direction = np.randon.rand()
#       position = np.random.randint(self.board_size, size=2)
#       rotquad90(self.board, direction, position)
#   if not self.step_count % self.steps_per_sample:
#       self.sample()
#   if not self.step_count % self.steps_per_board_sample:
#       board_sample_num = self.step_count // self.steps_per_board_sample
#       self.samples_board[board_sample_num] = self.board
    self.step_count += 1

## process data

def stratificatied(self):
    res = np.empty((self.samples_num, self.genotype_num))
    for i in range(self.genotype_num):
        res[:,i] = np.array([self.samples_frequency[:,i] + np.sum(self.samples_frequency[:,:i])])
    return res
        
def imagify_data(self):
    ## package boards' data into a displayable array.
    return np.array([self.S, self.R, self.C])

def display_frequency_timeseries(self):
    for i in range(8):
        pl.plot(np.arange(self.samples_num), self.samples_frequency[:,i], label=str(i), fillstyle='bottom')

property step_count:
    def __get__(self):
        return self.step_count
property generations:
    def __get__(self):
        return self.generations
    


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

cdef rotquad90(char[:,:,:] board, int direction, int[:] position):
    """
    rotquad90(self, direction, position) -> NoneType

    Turns a 2 by 2 sub-array of self.board by 90 degrees.
    Direction:
        0 - turn anticlockwise
        1 - turn clockwise
    with its lowest valued coordinate (upper left) in the "position" coordinate value.
    """

    cdef int temp_value
    cdef int row_i, col_i, gene_i
    row_i, col_i, gene_i = 0, 0, 0
    cdef int row0, col0
    row0, col0 = ( position[0],
                   position[1] )

    cdef int row1, col1
    row1, col1 = ( (row0 + 1) % board.shape[0],
                   (col0 + 1) % board.shape[1] )

    if direction == 0:
        for gene_i in range(3):
            ( board[row0, col0, gene_i],
              board[row0, col1, gene_i],
              board[row1, col0, gene_i],
              board[row1, col1, gene_i] ) = ( board[row0, col1, gene_i],
                                              board[row1, col1, gene_i],
                                              board[row0, col0, gene_i],
                                              board[row1, col0, gene_i] )
                                              # 00 01 -> 01 11
                                              # 10 11 -> 00 10
    else:
        for gene_i in range(3):
            ( board[row0, col0, gene_i],
              board[row0, col1, gene_i],
              board[row1, col0, gene_i],
              board[row1, col1, gene_i] ) = ( board[row1, col0, gene_i],
                                              board[row0, col0, gene_i],
                                              board[row1, col1, gene_i],
                                              board[row0, col1, gene_i] )
                                              # 00 01 -> 10 00
                                              # 10 11 -> 11 01
