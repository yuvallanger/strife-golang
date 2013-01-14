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
import pylab as pl
import timeit
import sys
import numba

usage = '''
{0}: A model

{0} <-h>
    <-d/--datefile> [data_filename.h5] <-c/--config> [config_filename]
'''

RECEPTOR = 0
SIGNAL = 1
COOPERATION = 2

def load_config(config_filename, default_config):
    """
    Takes a string holding filename and returns a dictionary with all the configuration values.
    """
    our_config = default_config
    if config_filename is None:
        return our_config

    config = ConfigParser.SafeConfigParser()
    config.read(config_filename)

    for key, val in our_config.items():
        our_config[key] = config.get('Config', key)

    return our_config

default_config = {
        'S_COST':                1,   # Metabolic cost of signalling
        'R_COST':                3,   # Metabolic cost of having a receptor
        'C_COST':                30,  # Metabolic cost of being cooperative
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


data_filename = default_config['data_filename']
conf_filename = None

default_config_filename = 'default.conf'
if os.path.exists(default_config_filename):
    conf_filename = default_config_filename

for i, arg in enumerate(sys.argv):
    if arg in ('--help', '-h'):
        raise usage
    if arg in ('--data', '-d'):
        print arg, sys.argv[i+1]
        data_filename = sys.argv[i+1]
    if arg in ('--config', '-c'):
        print arg, sys.argv[i+1]
        conf_filename = sys.argv[i+1]

config = load_config(conf_filename, default_config)
config['data_filename'] = data_filename

if os.path.exists(config['data_filename']):
    data = np.load(config['data_filename']).items()

    S_COST = data['S_COST']
    R_COST = data['R_COST']
    C_COST = data['C_COST']
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
    go()
else:
    #########
    # model parameters,
    #########

    #######
    # Cost of gene expression
    S_COST = int(config['S_COST'])
    R_COST = int(config['R_COST'])
    C_COST = int(config['C_COST'])

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

    ## Probability of each single diffusion operation
    diffusion_amount = float(config['diffusion_amount'])

    #######
    # radius of Signal or Cooperation effects.
    S_rad = int(config['S_rad'])
    G_rad = int(config['G_rad'])

    generations = int(config['generations'])

    print globals()['board_size']
    print int(config['board_size'])
    globals()['board_size'] = int(config['board_size'])
    print globals()['board_size']
    print int(config['board_size'])
    raw_input()

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

    # TODO do whole board sum
    # S_sum = np.empty((board_size, board_size),

    genotype_num = int(8)

    ## data sampling
    # we will take a frequency sample some number of times per generation
    steps_per_gen = int(board_size ** 2)
    samples_per_gen = int(1)
    samples_num = int(samples_per_gen * generations)
    samples_board_num = int(samples_num // 10)
    steps_per_sample = int(np.floor(1.0 * steps_per_gen // samples_per_gen))
    steps_per_board_sample = int(10 * steps_per_sample)
    sample_count = int(0)
    # We want to know the frequency of each genotype per generation
    samples_frequency = np.empty((samples_num, genotype_num), dtype=int)
    samples_nhood = np.empty((samples_num, genotype_num, genotype_num), dtype=int)
    samples_board = np.empty((samples_board_num, board_size, board_size, 3), dtype=np.int8)
    np.savez(config['data_filename'], **model.get_state())

###########################

S_COST = R_COST = C_COST = S_th = G_th = S_rad = G_rad = generations = \
genotype_num = steps_per_gen = \
samples_per_gen = samples_num = samples_board_num = steps_per_sample = \
steps_per_board_sample = sample_count = \
        board_size = step_count = steps_final = int(0)

baseline_cost = benefit = mutation_rate_r = mutation_rate_s = mutation_rate_c = \
        diffusion_amount = 1.0

board = np.empty((1,1,1), dtype='int8')
samples_frequency = np.empty((1,1), dtype='int32')
samples_nhood = np.empty((1,1,1), dtype='int32')
samples_board = np.empty((1,1,1,1), dtype='int8')

def init(config):
    global S_COST, R_COST, C_COST, baseline_cost, benefit, mutation_rate_r, \
           mutation_rate_s, mutation_rate_c, S_th, G_th, diffusion_amount, \
           S_rad, G_rad, generations, board_size, step_count, steps_final, \
           board, genotype_num, steps_per_gen, \
           samples_per_gen, samples_num, samples_board_num, steps_per_sample, \
           steps_per_board_sample, sample_count, samples_frequency, \
           samples_nhood, samples_board

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


@numba.jit('void(int8[:,:,:],int32,int32,float64,float64,float64)')
def mutate(board, pos_row, pos_col, p_r, p_s, p_c):
    raw_input((type(pos_row), type(pos_col), type(p_r), type(p_s), type(p_c),board_size))
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


def get_state():
    """
    The state of the simulation in a key, val dictionary.
    """
    return {'S_COST': S_COST,
            'R_COST': R_COST,
            'C_COST': C_COST,
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
    global S_COST, R_COST, C_COST, baseline_cost, benefit, mutation_rate_r, mutation_rate_s, mutation_rate_c, S_th, G_th, diffusion_amount, S_rad, G_rad, generations, board_size, step_count, steps_final, board, genotype_num, steps_per_gen, samples_per_gen, samples_num, samples_board_num, steps_per_sample, steps_per_board_sample, sample_count, samples_frequency, samples_nhood, samples_board

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



def rel_pos_find(i, axis):
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


#@numba.jit('i8[:,:,:](i8[:,:,:])')
def competitionstep(board):
    pos_1s = np.random.randint(board_size, size=(board_size ** 2, 2))
    rel_pos_s = np.random.randint(4, size=(board_size ** 2))
    p_pairs = np.random.rand(board_size ** 2, 2)
    p_muts = np.random.rand(board_size ** 2, 3)
    for i in range(board_size ** 2):
        print i
        pos_2_row = pos_1s[i, 0] + rel_pos_find(rel_pos_s[i], 0)
        pos_2_col = pos_1s[i, 1] + rel_pos_find(rel_pos_s[i], 1)
        pos_2t_row = pos_2_row % board_size
        pos_2t_col = pos_2_col % board_size

        if (board[pos_1s[i, 0], pos_1s[i, 1], :] == board[pos_2t_row, pos_2t_col, :]).all():
            board = mutate(board, pos_1s[i, 0],
                   pos_1s[i, 1],
                   p_muts[i, 0],
                   p_muts[i, 1],
                   p_muts[i, 2])
            continue

        winner = competition(pos_1s[i, 0],
                                  pos_1s[i, 1],
                                  pos_2_row,
                                  pos_2_col,
                                  p_pairs[i, 0],
                                  p_pairs[i, 1])
        if winner == 1:
            board[pos_2_row, pos_2_col, :] = board[pos_1s[i, 0], pos_1s[i, 1], :]
            board = mutate(board, pos_2_row,
                        pos_2_col,
                        p_muts[i, 0],
                        p_muts[i, 1],
                        p_muts[i, 2])
        else:
            board = copycell(board, pos_2_row,
                          pos_2_col,
                          pos_1s[i, 0],
                          pos_1s[i, 1])
            board = mutate(board, pos_1s[i, 0],
                        pos_1s[i, 1],
                        p_muts[i, 0],
                        p_muts[i, 1],
                        p_muts[i, 2])
    return board

#@numba.jit('i8[:,:,:](i8[:,:,:])')
def nextstep(board): #board, step_count):
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
    board = competitionstep(board)

    for i in range(int((board_size ** 2) * diffusion_amount) // 4):
        direction = np.random.randint(2)
        position = np.random.randint(board_size, size=2)
        board = rotquad90(board, direction, position)
    return board

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


@numba.jit('void(int8[:,:,:],int32,int32[:])')
def rotquad90(board, direction, position):
    """
    rotquad90(self, direction, position) -> NoneType

    Turns a 2 by 2 sub-array of self.board by 90 degrees.
    Direction:
        0 - turn anticlockwise
        1 - turn clockwise
    with its lowest valued coordinate (upper left) in the "position" coordinate value.
    """

    row0 = position[0]
    col0 = position[1]

    row1 = (row0 + 1) % board.shape[0]
    col1 = (col0 + 1) % board.shape[1]

    if direction == 0:
        for gene_i in range(3):
            board[row0, col0, gene_i] = board[row0, col1, gene_i]
            board[row0, col1, gene_i] = board[row1, col1, gene_i]
            board[row1, col0, gene_i] = board[row0, col0, gene_i]
            board[row1, col1, gene_i] = board[row1, col0, gene_i]
                                              # 00 01 -> 01 11
                                              # 10 11 -> 00 10
    else:
        for gene_i in range(3):
            board[row0, col0, gene_i] = board[row1, col0, gene_i]
            board[row0, col1, gene_i] = board[row0, col0, gene_i]
            board[row1, col0, gene_i] = board[row1, col1, gene_i]
            board[row1, col1, gene_i] = board[row0, col1, gene_i]
                                              # 00 01 -> 10 00
                                              # 10 11 -> 11 01

print 'go!' * 10
#    t = time.time()
#    every = 30*60
# TODO: Maybe add f and d somehow like in printf? {0}f {1}d
#    print "t: {0:f}, steps thus far: {1}".format(t, a.step_count)
#    steps_a = a.step_count
for step_count in range(5):  # range(model.generations):
    print step_count
    board = nextstep(board)
    if not step_count % steps_per_sample:
        sample()
    if not step_count % steps_per_board_sample:
        board_sample_num = step_count // steps_per_board_sample
        samples_board[board_sample_num] = board
    step_count += 1
#        delta_t = time.time() - t
#        if delta_t > every:
#            t = time.time()
#            a.save_h5()
#            steps_delta = a.step_count - steps_a
#            steps_a = a.step_count
#            eta = 1.0 * delta_t / (steps_delta+1) * (a.steps_final - a.step_count)
#            print "t: {0:f}, approx. time to fin: {1:f}".format(t, eta)
#            print "steps taken = {0}, steps since last save = {1}".format(a.step_count, steps_delta)
#            sys.exit(1)

np.savez(config['data_filename'], **model.get_state())
