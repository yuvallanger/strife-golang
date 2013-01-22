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

@numba.jit('i8[:](u8[:],i8,i8)')
def pymt64randint(mt, b, n):
    """
    pymt64randint(mt, b, n) -> [0,b), shape=(n,)
    """
    return np.int64(np.floor(pymt64.uniform(mt, n) * b))
    
def load_config(conf_filename, default_config):
    """
    Takes a string holding filename and returns a dictionary with all the configuration values.
    """
    our_config = default_config
    print 'load_config(),', our_config, conf_filename
    if conf_filename is None:
        return our_config

    config = ConfigParser.SafeConfigParser()
    config.read(conf_filename)

    for key, val in our_config.items():
        print 'before load_config(),', key, val
        our_config[key] = config.get('Config', key)
        print 'after load_config(),', key, val

    print 'end load_config(),', our_config

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
    print "os.path.exists(" + default_conf_filename + ") == True"
    conf_filename = default_conf_filename

for i, arg in enumerate(sys.argv):
    print 'sys.argv ==', i, arg
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
    print 'data_filename is None and conf_filename is None' 
    data_filename = default_config['data_filename']

def get_randomstate():
    return {'randomstate_current_' + str(i): val for i, val in enumerate(np.random.get_state())}

def get_state():
    """
    The state of the simulation in a key, val dictionary.
    """
    print 'get_state()'
    state =  {'S_cost': S_cost,
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
              'board_strain': board_strain,
              'board_signal_num': board_signal_num,
              'board_pg_num': board_pg_num,
              'board_prod': board_prod,
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
              'diffusion_step_num': diffusion_step_num}
    state.update({'randomstate_current_' + str(i): val for i, val in enumerate(np.random.get_state())}) 
    state.update({'randomstate_start_' + str(i): val for i, val in enumerate(np.random.get_state())})
    return state



if os.path.exists(config['data_filename']):
    print 'os.path.exists(' + config['data_filename'] + ')'
    data = {key: val for key, val in np.load(config['data_filename']).items()}

    #mt_0 = data['mt_0']
    #mt = data['mt']
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
    diffusion_step_num = data['diffusion_step_num']
    randomstate_start = (data['randomstate_start_0'],
                         data['randomstate_start_1'],
                         data['randomstate_start_2'],
                         data['randomstate_start_3'],
                         data['randomstate_start_4'])
    np.random.set_state(data['randomstate_current_0'],
                        data['randomstate_current_1'],
                        data['randomstate_current_2'],
                        data['randomstate_current_3'],
                        data['randomstate_current_4'])
    np.random.set_state(randomstate)
else:
    print 'os.path.exists(' + config['data_filename'] + ')'
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

    np.random.seed(0)
    randomstate_0 = np.random.get_state()

    # A cell can be Signalling and/or Receptive and/or Cooperative
    R = np.random.rand(board_size, board_size) < float(config['initial_receptives_amount'])
    S = np.random.rand(board_size, board_size) < float(config['initial_signallers_amount'])
    R = np.int8(R)
    S = np.int8(S)
    print R
    print S
    print C
    board_strain = R + 2 * S
    print "board_strain", type(board_strain), board_strain.dtype
    for i in range(3):
        print board[:,:,i]

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


    randomrng = random.Random(0)

    np.savez(config['data_filename'], **get_state())

###########################


######
## functions
######

def signal_count(pos_row, pos_col, signal):
    signal_sum = 0
    for s_row_i in range(pos_row - S_rad, pos_row + S_rad + 1):
        for s_col_i in range(pos_col - S_rad, pos_col + S_rad + 1):
            if board[s_row_i % board_size, s_col_i % board_size, SIGNAL] == signal:
                signal_sum += 1
    return signal_sum


#@numba.autojit
@numba.jit('i8(uint8[:,:,:],i8,i8)')
def goods_count(board, pos_row, pos_col):
    """
    Counts the number of public goods at a certain position
    """
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
            #print (g_row_i, g_col_i), board[g_row_i%board_size, g_col_i%board_size],
            #print goods_sum, 
            # Able to produce public goods
            if int(board[g_row_i % board_size,
                     g_col_i % board_size,
                     COOPERATION]) == 1:
                #print 'C,',
                # Isn't receptive
                if int(board[g_row_i % board_size,
                         g_col_i % board_size,
                         RECEPTOR]) == 0:
                    #print 'R,',
                    goods_sum += 1
                # Receptive and signal reaches signal threshold
                elif int(signal_count(g_row_i, g_col_i, 1)) >= S_th:
                    #print 'S_th+',
                    goods_sum += 1
                #else:
                    #print 'S_th-',
            #else:
                #print 'c',
            #print goods_sum
    return goods_sum



#@numba.autojit
@numba.jit('f8(i8,i8)')
def fitness(pos_row, pos_col):
    pos_row_t = pos_row % board_size
    pos_col_t = pos_col % board_size

    #print 'fitness(),', 'pos_row_t, ', pos_row_t
    #print 'fitness(),', 'pos_col_t, ', pos_col_t

    goods_num = goods_count(board, pos_row_t, pos_col_t)
    #print S_th, G_th
    #board_show_r = np.array([[board.transpose(2,0,1)[0,row_i, col_i] for col_i in np.arange(pos_col-2, pos_col+3) % board_size] for row_i in np.arange(pos_row-2, pos_row+3) % board_size])
    #board_show_s = np.array([[board.transpose(2,0,1)[1,row_i, col_i] for col_i in np.arange(pos_col-2, pos_col+3) % board_size] for row_i in np.arange(pos_row-2, pos_row+3) % board_size])
    #board_show_c = np.array([[board.transpose(2,0,1)[2,row_i, col_i] for col_i in np.arange(pos_col-2, pos_col+3) % board_size] for row_i in np.arange(pos_row-2, pos_row+3) % board_size])
    #print 'fitness(),', 'goods_num,', goods_num
    #print 'fitness(),', 'receptor:'
    #print board_show_r
    #print 'fitness(),', 'signal:'
    #print board_show_s
    #print 'fitness(),', 'cooperation:'
    #print board_show_c
    #raw_input()

    result = benefit * (goods_num >= G_th)
    result = 1 - result
    result = (baseline_cost + result) / result

    #print 'fitness(),', 'result,', result#, stuffprint(result)
    return result

def stuffprint(a):
    print 'stuffprint(),', type(a), a.dtype, a.shape


#@numba.jit('i8[:,:](uint8[:,:,:], i8,i8,i8,i8, f8,f8)')
@numba.jit('i8(uint8[:,:,:],i8,i8,i8,i8,f8,f8)')
def competition(board,
                cell_a_row,
                cell_a_col,
                cell_b_row,
                cell_b_col,
                cell_a_p,
                cell_b_p):
    """
    competition(pair_i, pos_1s, pos_2s, p_s) -> (int, int)

    Decides which of the two positions wins.

    Coordinates are a numpy array of shape=(2,) and an integer dtype.
    Takes two adjacent position coordinates on the board and each one's TODO: what's the name of such probability values?
    p_pair: an array of two uniform distribution over [0, 1).
    """

    score_a = cell_a_p * fitness(cell_a_row, cell_a_col)

    score_b = cell_b_p * fitness(cell_b_row, cell_b_col)

    if (float(score_a) > float(score_b)):
        return 0
        #return (cell_a_row, cell_a_col, cell_b_row, cell_b_col)
    else:
        return 1
        #return (cell_a_row, cell_a_col, cell_b_row, cell_b_col)



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


#@numba.jit('void(uint8[:,:,:],i8,i8,f8,f8,f8)')
def mutate(board, pos_row, pos_col, p_r, p_s, p_c):
    """
    mutate(board, self, pos) -> NoneType

    For each value of self.board at position "pos",
        change its value at probability of self.mutation_rate_[r/s/c].
    """
    pos_row_t = pos_row % board_size
    pos_col_t = pos_col % board_size
    if float(p_r) < float(mutation_rate_r):
        if board[pos_row_t, pos_col_t, 0]:
            board[pos_row_t, pos_col_t, 0] = 0
        else:
            board[pos_row_t, pos_col_t, 0] = 1
    if float(p_s) < float(mutation_rate_s):
        if board[pos_row_t, pos_col_t, 1]:
            board[pos_row_t, pos_col_t, 1] = 0
        else:
            board[pos_row_t, pos_col_t, 1] = 1
    if float(p_c) < float(mutation_rate_c):
        if board[pos_row_t, pos_col_t, 2]:
            board[pos_row_t, pos_col_t, 2] = 0
        else:
            board[pos_row_t, pos_col_t, 2] = 1
    #return board


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


#@numba.jit('i8(i8[:,:],i8[:,:])')
#@numba.autojit
def rel_pos_find(a, directions):
    """
    Takes two (n, 2) integer arrays (*a* and *b*), and according to a (n,) integers array holding the numbers [0,4] denoting directions.
    """
    d = np.array([[-1,0], [0, -1], [1, 0], [0, 1]], dtype=np.int64) 
    return a + d[directions]
    #for i in range(a.shape[0]):
    #    if directions[i] == 0:
    #        print 'a'
    #    else:
    #        print 'b'
    #        b[i, 0] = (a[i, 0] - 1) % board_size
    #        b[i, 1] = a[i, 1]
    #    elif (directions[i] == 1):
    #        b[i, 0] = a[i, 0]
    #        b[i, 1] = (a[i, 1] - 1) % board_size
    #    elif (directions[i] == 2):
    #        b[i, 0] = (a[i, 0] + 1) % board_size
    #        b[i, 1] = a[i, 1]
    #    else:
    #        b[i, 0] = a[i, 0]
    #        b[i, 1] = (a[i, 1] + 1) % board_size
    #code = '''
    #int i = 0;
    #for (i = 0; i < Na[0]; i++) {
    #    if (DIRECTIONS1(i) == 0) {
    #        B2(i, 0) = pyModulus(A2(i, 0) - 1, board_size);
    #        B2(i, 1) = A2(i, 1);
    #        continue;
    #    }
    #    if (DIRECTIONS1(i) == 1) {
    #        B2(i, 0) = A2(i, 0);
    #        B2(i, 1) = pyModulus(A2(i, 1) - 1, board_size);
    #        continue;
    #    }
    #    if (DIRECTIONS1(i) == 2) {
    #        B2(i, 0) = pyModulus(A2(i, 0) + 1, board_size);
    #        B2(i, 1) = A2(i, 1);
    #        continue;
    #    }
    #    B2(i, 0) = A2(i, 0);
    #    B2(i, 1) = pyModulus(A2(i, 1) + 1, board_size);
    #}
    #'''

    #support_code = r'''
    #int pyModulus(int a, int b) {
    #    return ((a % b) + b) % b;
    #}
    #'''

    #sp.weave.inline(code, ['a', 'b', 'directions', 'board_size'], support_code=support_code)

#@numba.autojit
def same(board, pair_i, pos_1s, pos_2s):
    for gene_i in range(3):
        if board[pos_1s[pair_i, 0], pos_1s[pair_i, 1], gene_i] != board[pos_2s[pair_i, 0], pos_2s[pair_i, 1], gene_i]:
            return 0
    return 1

#@numba.jit('void(uint8[:,:,:])')
def competitionstep(board):
    pos_1s = np.random.randint(board_size, size=(board_size ** 2, 2))
    rel_pos_s = np.random.randint(4, size=(board_size ** 2))
    p_pairs = np.random.rand(board_size ** 2, 2)
    p_muts = np.random.rand(board_size ** 2, 3)
    pos_2s = rel_pos_find(pos_1s, rel_pos_s)
    pos_2st = pos_2s % board_size
    a_rs, a_cs, b_rs, b_cs = pos_1s[:, 0], pos_1s[:, 1], pos_2s[:, 0], pos_2s[:, 1]
    
    for pair_i in range(board_size ** 2):
        #print 'competitionstep(), pair_i first,', pair_i
        equal = sp.weave.inline(r'''
                           int gene_i;
                           return_val = 1;
                           for (gene_i = 0; gene_i < 3; gene_i++) {
                               if (BOARD3(a_r, a_c, gene_i) != BOARD3(b_r, b_c, gene_i)) {
                                   return_val = 0;
                                   break;
                               }
                           }''',
                  ['a_r', 'a_c', 'b_r', 'b_c', 'board'],
                  {'a_r': int(pos_1s[pair_i, 0]),
                   'a_c': int(pos_1s[pair_i, 1]),
                   'b_r': int(pos_2st[pair_i, 0]),
                   'b_c': int(pos_2st[pair_i, 1]),
                   'board': board})
        #print 'competitionstep(), pair_i equal,', pair_i
        #print 'competitionstep(),', equal, board[pos_1s[pair_i], pos_1s[pair_i], :], board[pos_2s[pair_i], pos_2s[pair_i], :]
        if equal:
            #print 'competitionstep(), pair_i if equal, ', pair_i
            mutate(board,
                   pos_1s[pair_i, 0],
                   pos_1s[pair_i, 1],
                   p_muts[pair_i, 0],
                   p_muts[pair_i, 1],
                   p_muts[pair_i, 2])
            #print 'competitionstep(), pair_i if equal end, ', pair_i
        else:
            #print 'competitionstep(), pair_i if not equal, ', pair_i
            competition(board,
                                                                         pos_1s[pair_i, 0],
                                                                         pos_1s[pair_i, 1],
                                                                         pos_2s[pair_i, 0],
                                                                         pos_2s[pair_i, 1],
                                                                         p_pairs[pair_i, 0],
                                                                         p_pairs[pair_i, 1])
            #print 'competitionstep(), pair_i if not equal end, ', pair_i
    #sp.weave.inline(r'''
    #for (pair_i = 0; pair_i < cells_num; pair_i)
    #{
    #    int ret_val;
    #    int gene_i = 0;
    #    //printf("pair_i: %d, ", pair_i);
    #    return_val = ret_val;
    #    for (gene_i = 0; gene_i < 3; gene_i++)
    #    {
    #        // printf("gene_i: %d, ", gene_i);
    #        if (BOARD3(POS_1S2(pair_i, 0), POS_1S2(pair_i, 1), gene_i) != BOARD3(POS_2S2(pair_i, 0), POS_2S2(pair_i, 1), gene_i))
    #        {
    #            return_val = 0;
    #            //printf("!=, ");
    #            break;
    #        }
    #    }
    #    return_val = ret_val;
    #    //printf("ret_val: %d ", ret_val);
    #    //printf("\n");
    #''', ['pair_i', 'board', 'pos_1s', 'pos_2s'])
    #    #if same(board, pair_i, pos_1s, pos_2s):
    #    #    pass
    #    #    #mutate(board, pos_1s[i, 0],
    #    #    #       pos_1s[i, 1],
    #    #    #       p_muts[i, 0],
    #    #    #       p_muts[i, 1],
    #    #    #       p_muts[i, 2])
    #    #    #continue
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
def rotquad90(board, directions, positions):
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


@numba.jit('void(uint8[:,:,:])') # , uint64[:])')
def nextstep(board):
    competitionstep(board)
    directions = np.random.randint(2, size=diffusion_step_num)
    positions = np.random.randint(board_size, size=(diffusion_step_num, 2))
    rotquad90(board, directions, positions)

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


@numba.jit('void(int8[:,:,:],f8[:])')
def go(board, times):
    print 'start go()'
    #    every = 30*60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    #    steps_a = a.step_count
    print 'go(),', board_size, generations
    times[0] = np.float64(time.time())
    for step_count in range(generations):
        nextstep(board)
        times[step_count + 1] = time.time()
        #if (step_count % 25) == 0:
        print 'go(),', step_count
        print 'go(),', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 300**2 / 60 / 60
        print 'go(),', np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 300**2 / 60 / 60
        print 'go(),', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 150**2 / 60 / 60
        print 'go(),', np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 150**2 / 60 / 60
        print 'go(),', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 75**2 / 60 / 60
        print 'go(),', np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 75**2 / 60 / 60
        print 'go(),', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 50**2 / 60 / 60
        print 'go(),', np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 50**2 / 60 / 60
        print 'go(),', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 25**2 / 60 / 60
        print 'go(),', np.average(times[1:step_count + 2] - times[:step_count+1]) * 10000 / board_size**2 * 25**2 / 60 / 60
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
        

times = np.empty((generations+1), dtype=np.float64)

print 'root level,', times.dtype, times.shape, type(times)
print 'root level, before go()'
go(board, times)

l=[]
for time_0, time_1 in zip(times[:-1], times[1:]):
    l.append(time_1-time_0)
print 'root level,', sum(l) / generations * 10000 / 10**2 * 300**2 / 60 / 60 

    
np.savez(config['data_filename'], **get_state())
