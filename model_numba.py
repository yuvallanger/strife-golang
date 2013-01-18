# -*- coding: utf-8 -*-
# cython: profile=True

"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

from __future__ import division
import os
import time
import scipy as sp
#cimport scipy as sp
import scipy.signal
import scipy.weave
import numpy as np
import pylab as pl
import sys
import numba
import ConfigParser
import pymt64
import random

usage = '''
{0}: A model

{0} <-h>
    <-d/--datefile> [data_filename.h5] <-c/--config> [conf_filename]
'''

RECEPTOR = 0
SIGNAL = 1
COOPERATION = 2

def load_config(conf_filename, default_config):
    """
    Takes a string holding filename and returns a dictionary with all the configuration values.
    """
    our_config = default_config
    print our_config, conf_filename
    if conf_filename is None:
        return our_config

    config = ConfigParser.SafeConfigParser()
    config.read(conf_filename)

    for key, val in our_config.items():
        print key, val
        our_config[key] = config.get('Config', key)
    print our_config

    return our_config

default_config = {
        'S_cost':                1,   # Metabolic cost of signalling
        'R_cost':                3,   # Metabolic cost of having a receptor
        'C_cost':                30,  # Metabolic cost of being cooperative
        'baseline_cost':         100, # Basal metabolic cost

        'benefit':               0.9, # The fraction reduced, out of total metabolic cost, when
                                      #    public goods reach threshold of benefit.

        'mutation_rate_r': 1e-4, # Likelihoods of switch (on to off and vica versa) for each gene per cell per generation.
        'mutation_rate_s': 1e-4, #
        'mutation_rate_c': 1e-4, #

        'S_th':               3, # Amount of signal needed for a receptive and cooperative cell to
                                 #    start producing public goods.
        'G_th':               3, # Amount of public goods needed for metabolic benefit.

        'diffusion_amount': 0.5, # Average fraction of cells out of total cells on the board (board_size**2)
                                 #    which will be involved
        'board_size':        50, # The length of the board. The board will be of shape (board_size, board_size).
        'generations':       10, # Each generation involves, on average, all cells in a competition, meaning
                                 #    board_size**2 competitions.
        'S_rad':              1, # Radius of the signal's effect.
        'G_rad':              1,
        'samples_per_gen':    1,
        'initial_receptives_amount': 0.5,
        'initial_signallers_amount': 0.5,
        'initial_cooperators_amount': 0.5,
        'data_filename': 'default.npz'}


conf_filename = None
data_filename = None
default_conf_filename = 'default.conf'
if os.path.exists(default_conf_filename):
    print "default_conf_filename exists", default_conf_filename
    conf_filename = default_conf_filename

for i, arg in enumerate(sys.argv):
    print i, arg
    if arg in ('--help', '-h'):
        print "arg in ('--help', '-h')", arg
        raise usage
    if arg in ('--data', '-d'):
        print "arg in ('--data', '-d')", arg, sys.argv[i+1]
        data_filename = sys.argv[i+1]
    if arg in ('--config', '-c'):
        print "arg in ('--config', '-c')", arg, sys.argv[i+1]
        conf_filename = sys.argv[i+1]

config = load_config(conf_filename, default_config)
if data_filename is None and conf_filename is None:
        data_filename = default_config['data_filename']

def get_state():
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
            'samples_board': samples_board,
            'mt': mt}


if os.path.exists(config['data_filename']):
    print "there's a data file", config['data_filename']
    data = {key: val for key, val in np.load(config['data_filename']).items()}

    mt = data['mt']
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
else:
    print "there isn't a data file", config['data_filename']
    #########
    # model parameters,
    #########

    #######
    # Cost of gene expression
    S_cost = int(config['S_cost'])
    R_cost = int(config['R_cost'])
    C_cost = int(config['C_cost'])

    # B for Baseline, basal, "basic metabolic burden"
    baseline_cost = int(config['baseline_cost'])

    #####
    # benefit from cooperation. "reward factor" in the article.
    benefit = float(config['benefit'])

    ######
    # mutation per generation
    mutation_rate_r = float(config['mutation_rate_r'])
    mutation_rate_s = float(config['mutation_rate_s'])
    mutation_rate_c = float(config['mutation_rate_c'])

    ## neighbours effects' thresholds
    S_th = int(config['S_th'])
    # quorum threshold
    G_th = int(config['G_th'])
    # Cooperation threshold. Above it, public goods makes a difference.
    
    board_size = int(config['board_size'])

    ## Probability of each single diffusion operation
    diffusion_amount = float(config['diffusion_amount'])
    diffusion_step_num = int((board_size ** 2) * diffusion_amount) // 4

    #######
    # radius of Signal or Cooperation effects.
    S_rad = int(config['S_rad'])
    G_rad = int(config['G_rad'])

    generations = int(config['generations'])


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
    step_count = int(0)


    steps_final = int(generations * board_size**2)

    # A cell can be Signalling and/or Receptive and/or Cooperative
    R = np.random.rand(board_size, board_size) < float(config['initial_receptives_amount'])
    S = np.random.rand(board_size, board_size) < float(config['initial_signallers_amount'])
    C = np.random.rand(board_size, board_size) < float(config['initial_cooperators_amount'])
    R = np.int8(R)
    S = np.int8(S)
    C = np.int8(C)
    board = np.array([R, S, C]).transpose((1,2,0))
    print "board", type(board), board.dtype

    # TODO do whole board sum
    # S_sum = np.empty((board_size, board_size),

    genotype_num = int(8)

    ## data sampling
    # we will take a frequency sample some number of times per generation
    steps_per_gen = int(board_size ** 2)
    samples_per_gen = int(1)
    samples_num = int(samples_per_gen * generations)
    samples_board_num = int(samples_num // 100)
    steps_per_sample = int(np.floor(1.0 * steps_per_gen // samples_per_gen))
    steps_per_board_sample = int(10 * steps_per_sample)
    sample_count = int(0)
    # We want to know the frequency of each genotype per generation
    samples_frequency = np.empty((samples_num, genotype_num), dtype=int)
    samples_nhood = np.empty((samples_num, genotype_num, genotype_num), dtype=int)
    samples_board = np.empty((samples_board_num, board_size, board_size, 3), dtype=np.int8)

    mt = pymt64.init(0)
    print mt.dtype

    randomrng = random.Random(0)

    np.savez(config['data_filename'], **get_state())

###########################


######
## functions
######

def signal_count(row, col, signal):
    signal_sum = 0
    for s_row_i in range(row - S_rad, row + S_rad + 1):
        for s_col_i in range(col - S_rad, col + S_rad + 1):
            if board[s_row_i % board_size, s_col_i % board_size, SIGNAL] == signal:
                signal_sum += 1
    return signal_sum


def goods_count(pos_row, pos_col):
    """
    Counts the number of public goods at a certain position
    """
    #cdef int s_row_i, s_col_i, g_row_i, g_col_i, goods_sum

#     has fitness effect?:
#         surrounded by enough goods producing cells
#     produces goods?:
#         no receptor and cooperator
#         or
#         surrounded by enough signal producing cells and
#             is receptive and
#             is cooperator
    goods_sum = 0
    for g_row_i in range(pos_row - G_rad,
                         pos_row + G_rad + 1):
        for g_col_i in range(pos_col - G_rad,
                             pos_col + G_rad + 1):
            # Able to produce public goods
            if board[g_row_i % board_size,
                     g_col_i % board_size,
                     COOPERATION] == 1:
                # Isn't receptive
                if board[g_row_i % board_size,
                         g_col_i % board_size,
                         RECEPTOR] == 0:
                    goods_sum += 1
                # Receptive and signal reaches signal threshold
                elif signal_count(g_row_i, g_col_i, 1) >= S_th:
                    goods_sum += 1
    return goods_sum



def fitness(pos_row, pos_col, p):
    pos_row_t, pos_col_t = (pos_row % board_size,
                            pos_row % board_size)

    result = benefit * (goods_count(pos_row_t, pos_col_t) >= G_th)
    result = 1 - result
    result = (baseline_cost + result) / result
    result = result * p

    return result


#numba.jit('void(i8,i8,i8,i8,float64,float64'))
def competition(c_pos_1_row,
                c_pos_1_col,
                c_pos_2_row,
                c_pos_2_col,
                p_1,
                p_2):
    """
    competition(self, cell_pos_1, cell_pos_2, p_pair) -> (int, int)

    Decides which of the two positions wins.

    Coordinates are a numpy array of shape=(2,) and an integer dtype.
    Takes two adjacent position coordinates on the board and each one's TODO: what's the name of such probability values?
    p_pair: an array of two uniform distribution over [0, 1).
    """
#    cdef double score_1, score_2
#    cdef int c_pos_2t_row, c_pos_2t_col
    c_pos_2t_row = c_pos_2_row % board_size
    c_pos_2t_col = c_pos_2_col % board_size

    score_1 = p_1 * fitness(c_pos_1_row,
                            c_pos_1_col,
                            p_1)

    score_2 = p_2 * fitness(c_pos_2t_row,
                            c_pos_2t_col,
                            p_2)

    return (1 if score_1 > score_2 else 2)



def copycell(board, orig_row, orig_col,
             dest_row, dest_col):
    """
    copycell(board, orig_row, orig_col, dest_row, dest_col) -> NoneType

    Copies the contents of self.board at coordinates of position "orig" into the position of coordinates "dest".
    Coordinates are a numpy array of shape=(2,) and an integer dtype.
    """
    #cdef int i, orig_row_t, orig_col_t, dest_row_t, dest_col_t
    orig_row_t = orig_row % board_size
    orig_col_t = orig_col % board_size
    dest_row_t = dest_row % board_size
    dest_col_t = dest_col % board_size

    for i in range(3):
        board[dest_row_t, dest_col_t, i] = board[orig_row_t, orig_col_t, i]
    return board


#@numba.jit('int8[:,:,:](int8[:,:,:],int32,int32,float64,float64,float64)')
def mutate(board, pos_row, pos_col, p_r, p_s, p_c):
    """
    mutate(board, self, pos) -> NoneType

    For each value of self.board at position "pos",
        change its value at probability of self.mutation_rate_[r/s/c].
    """
    pos_row_t = pos_row % board_size
    pos_col_t = pos_col % board_size
    if p_r < mutation_rate_r:
        if board[pos_row_t, pos_col_t, 0]:
            board[pos_row_t, pos_col_t, 0] = 0
        else:
            board[pos_row_t, pos_col_t, 0] = 1
    if p_s < mutation_rate_s:
        if board[pos_row_t, pos_col_t, 1]:
            board[pos_row_t, pos_col_t, 1] = 0
        else:
            board[pos_row_t, pos_col_t, 1] = 1
    if p_c < mutation_rate_c:
        if board[pos_row_t, pos_col_t, 2]:
            board[pos_row_t, pos_col_t, 2] = 0
        else:
            board[pos_row_t, pos_col_t, 2] = 1
    return board


def sample():
    global samples_nhood, samples_frequency, sample_count
    joint_board = board[:, :, 0] + 2 * board[:, :, 1] + 4 * board[:, :, 2]
    for gene in range(8):
        gene_board = joint_board == gene
        gene_frequency = np.sum(gene_board)
        # neighbours_genotype
        for nh_gene in range(8):
            nh_board = joint_board == nh_gene
            nh_gene_count = sp.signal.convolve2d(nh_board, np.ones((3,3)),
                                                 mode='same',
                                                 boundary='wrap')
            nh_gene_count_of_gene = np.sum(nh_gene_count * gene_board,
                                           dtype=int)
            samples_nhood[sample_count, gene, nh_gene] = nh_gene_count_of_gene
        samples_frequency[sample_count, gene] = gene_frequency
    sample_count += 1


#@numba.autojit
def rel_pos_find(a, b, directions):
    """
    Takes two (n, 2) integer arrays (*a* and *b*), and according to a (n,) integers array holding the numbers [0,4] denoting directions.
    """
    code = '''
    int i = 0;
    for (i = 0; i < Na[0]; i++) {
        if (DIRECTIONS1(i) == 0) {
            B2(i, 0) = pyModulus(A2(i, 0) - 1, board_size);
            B2(i, 1) = A2(i, 1);
            continue;
        }
        if (DIRECTIONS1(i) == 1) {
            B2(i, 0) = A2(i, 0);
            B2(i, 1) = pyModulus(A2(i, 1) - 1, board_size);
            continue;
        }
        if (DIRECTIONS1(i) == 2) {
            B2(i, 0) = pyModulus(A2(i, 0) + 1, board_size);
            B2(i, 1) = A2(i, 1);
            continue;
        }
        B2(i, 0) = A2(i, 0);
        B2(i, 1) = pyModulus(A2(i, 1) + 1, board_size);
    }
    '''

    support_code = r'''
    int pyModulus(int a, int b) {
        return ((a % b) + b) % b;
    }
    '''

    sp.weave.inline(code, ['a', 'b', 'directions', 'board_size'], support_code=support_code)

#@numba.jit('int64(uint64[:],int32,int32)')
def pymt64randint(mt, b, n):
    return np.int64(np.floor(pymt64.uniform(mt, n) * b))
    

@numba.autojit#('void(int8[:,:,:],uint64[:])')
def competitionstep(board, mt):
    pos_1s = pymt64randint(mt, board_size, 2 * board_size ** 2).reshape((board_size ** 2, 2))
    rel_pos_s = pymt64randint(mt, 4, board_size ** 2)
    p_pairs = pymt64.uniform(mt, 2 * board_size ** 2).reshape((board_size ** 2, 2))
    p_muts = pymt64.uniform(mt, 3 * board_size ** 2).reshape((board_size ** 2, 3))
    pos_2s = np.empty((board_size ** 2, 2), dtype=pos_1s.dtype)
    rel_pos_find(pos_1s, pos_2s, rel_pos_s)
    for i in range(board_size ** 2):
        pass
        #pos_1_row = np.random.randint(board_size) #, size=(board_size ** 2, 2))
        #pos_1_col = np.random.randint(board_size)
        #rel_pos = np.random.randint(4) #, size=(board_size ** 2))
        #pos_2t_row = pos_2_row % board_size
        #pos_2t_col = pos_2_col % board_size

        #if (board[pos_1s[i, 0], pos_1s[i, 1]] == board[pos_2t_row, pos_2t_col]).all():
        #    mutate(board, pos_1s[i, 0],
        #           pos_1s[i, 1],
        #           p_muts[i, 0],
        #           p_muts[i, 1],
        #           p_muts[i, 2])
        #    continue

        #winner = competition(pos_1s[i, 0],
        #                          pos_1s[i, 1],
        #                          pos_2_row,
        #                          pos_2_col,
        #                          p_pairs[i, 0],
        #                          p_pairs[i, 1])
    #    p_pair = np.random.rand(2)
    #    p_muts = np.random.rand(3)
#        if winner == 1:
#            board[pos_2_row, pos_2_col, :] = board[pos_1s[i, 0], pos_1s[i, 1], :]
#            mutate(board, pos_2_row,
#                        pos_2_col,
#                        p_muts[i, 0],
#                        p_muts[i, 1],
#                        p_muts[i, 2])
#        else:
#            copycell(board, pos_2_row,
#                          pos_2_col,
#                          pos_1s[i, 0],
#                          pos_1s[i, 1])
#            mutate(board, pos_1s[i, 0],
#                        pos_1s[i, 1],
#                        p_muts[i, 0],
#                        p_muts[i, 1],
#                        p_muts[i, 2])


#@numba.autojit#('void(int8[:,:,:],int32,int32[:])')
def rotquad90(board, diffusion_step_num, directions, positions):
    """
    rotquad90(self, direction, position) -> NoneType

    Turns a 2 by 2 sub-array of self.board by 90 degrees.
    Direction:
        0 - turn anticlockwise
        1 - turn clockwise
    with its lowest valued coordinate (upper left) in the "position" coordinate value.
    """
    code = r'''
    int row0, col0, row1, col1, i, gene_i, tmp; 
    for (i = 0; i < (int) diffusion_step_num; i++) {
        row0 = POS2(i, 0);
        col0 = POS2(i, 1);
        row1 = (row0 + 1) % Nboard[0];
        col1 = (col0 + 1) % Nboard[1];
        if DIRECTIONS1(i) {
            for (gene_i = 0; gene_i < 3; gene_i++) {
                tmp = BOARD3(row0, col0, gene_i);
                BOARD3(row0, col0, gene_i) = BOARD3(row0, col1, gene_i);
                BOARD3(row0, col1, gene_i) = BOARD3(row1, col1, gene_i);
                BOARD3(row1, col1, gene_i) = BOARD3(row1, col0, gene_i);
                BOARD3(row1, col0, gene_i) = tmp;
                                                  // # 00 01 -> 01 11
                                                  // # 10 11 -> 00 10
            }
        } else {
            for (gene_i = 0; gene_i < 3; gene_i++) {
                tmp = BOARD3(row0, col0, gene_i);
                BOARD3(row0, col0, gene_i) = BOARD3(row1, col0, gene_i);
                BOARD3(row1, col0, gene_i) = BOARD3(row1, col1, gene_i);
                BOARD3(row1, col1, gene_i) = BOARD3(row0, col1, gene_i);
                BOARD3(row0, col1, gene_i) = tmp;
            }
        }
    }
    '''

    sp.weave.inline(code,
                    ['pos', 'board', 'directions', 'diffusion_step_num'],
                    {'pos': positions,
                     'board': board,
                     'directions': directions,
                     'diffusion_step_num': diffusion_step_num})

    #if direction < 0.5:
    #    board[row0, col0, :], board[row0, col1, :], board[row1, col0, :], board[row1, col1, :] = board[row0, col1, :], board[row1, col1, :], board[row0, col0, :], board[row1, col0, :]
    #                                          # 00 01 -> 01 11
    #                                          # 10 11 -> 00 10
    #else:
    #    board[row0, col0, :], board[row0, col1, :], board[row1, col0, :], board[row1, col1, :] = board[row1, col0, :], board[row0, col0, :], board[row1, col1, :], board[row0, col1, :]
    #                                          # 00 01 -> 10 00
    #                                          # 10 11 -> 11 01


@numba.autojit#('void(int8[:,:,:],uint64[:])') # , uint64[:])')
def nextstep(board, mt):
#    global step_count, board
#    cdef int rel_pos_i, pos_1_row, pos_1_col, pos_2_row, pos_2_col, winner, pos_2t_row, pos_2t_col
#    cdef double p_1, p_2, p_r, p_s, p_c
#    cdef np.ndarray[np.int_t, ndim=2] pos_1s
#    cdef int[:,:] pos_1s_memview
#    cdef np.ndarray[np.int_t, ndim=1] rel_pos_s
#    cdef int[:] rel_pos_s_memview
#    cdef np.ndarray[np.double_t, ndim=2] p_pairs
#    cdef double[:,:] p_pairs_memview
#    cdef np.ndarray[np.double_t, ndim=2] p_muts
#    cdef double[:,:] p_muts_memview
#   pass 
    competitionstep(board, mt)
    directions = np.int8(np.floor(pymt64.uniform(mt, diffusion_step_num) * 2))
    positions = np.int16(np.floor(pymt64.uniform(mt, 2 * diffusion_step_num) * board_size)).reshape((diffusion_step_num, 2))
    rotquad90(board, diffusion_step_num, directions, positions)

## process data


#def stratificatied():
#    res = np.empty((self.samples_num, self.genotype_num))
#    for i in range(self.genotype_num):
#        res[:, i] = np.array([samples_frequency[:, i] +
#                             np.sum(samples_frequency[:, :i])])
#    return res
#
#
#def display_frequency_timeseries():
#    for i in range(8):
#        pl.plot(np.arange(samples_num),
#                samples_frequency[:, i],
#                label=str(i),
#                fillstyle='bottom')


#@numba.autojit#('void(i8[:,:,:],uint64[:],float64[:])')
def go(board, mt, times):
    #    every = 30*60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    #    steps_a = a.step_count
    print board_size, generations
    times[0] = time.time()
    for step_count in range(generations):
        print step_count
        nextstep(board, mt)
        times[step_count + 1] = time.time()
        print (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 300**2 / 60 / 60
        print np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 300**2 / 60 / 60
        #time.time()#1.0 * (time.time() - iter_start) / board_size**2 * 300**2 / 60 / 60 * 10000
#        if not step_count % steps_per_sample:
#            sample()
#        if not step_count % steps_per_board_sample:
#            board_sample_num = step_count // steps_per_board_sample
#            samples_board[board_sample_num] = board
#            delta_t = time.time() - t
#            if delta_t > every:
#                t = time.time()
#                a.save_h5()
#                steps_delta = a.step_count - steps_a
#                steps_a = a.step_count
#                eta = 1.0 * delta_t / (steps_delta+1) * (a.steps_final - a.step_count)
#                print "t: {0:f}, approx. time to fin: {1:f}".format(t, eta)
#                print "steps taken = {0}, steps since last save = {1}".format(a.step_count, steps_delta)
    print board
        

times = np.empty((generations+1), dtype=np.float64)


go(board, mt, times)

l=[]
for time_0, time_1 in zip(times[:-1], times[1:]):
    l.append(time_1-time_0)
print sum(l) / generations * 10000 / 10**2 * 300**2 / 60 / 60 

    
np.savez(config['data_filename'], **get_state())
