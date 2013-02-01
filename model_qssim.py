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

s4strain = np.array([0, 0, 1, 1], dtype='int8')  # 0: S0, 1: S0, 2: S1, 3: S1 
r4strain = np.array([0, 1, 0, 1], dtype='int8')  # 0: R0, 1: R1, 2: R0, 3: R1
strain_spec = np.array([[0, 2],
                        [1, 3]], dtype='int8')  # [r, s]
same_strain = np.array([[1 if i == j else 0 for j in range(4)] for i in range(4)])

#@numba.autojit
@profile
def init_board_signal_num(board_strain,
                          board_signal_num,
                          S_rad):
    support_code = r'''
    int pyModulus(int a, int b) {
        return ((a % b) + b) % b;
    }
    '''
    code = r'''
    int row_center_i, col_center_i, row_rad_i, col_rad_i, row_rad_i_t, col_rad_i_t;
    int s4strain[4] = {0, 0, 1, 1};  // # 0: S0, 1: S0, 2: S1, 3: S1 
    int r4strain[4] = {0, 1, 0, 1};  //, dtype='int8')  # 0: R0, 1: R1, 2: R0, 3: R1
    int center_signal_type;

    for (row_center_i = 0;
         row_center_i < Nboard_strain[0];
         row_center_i++)
    {
        for (col_center_i = 0;
             col_center_i < Nboard_strain[1];
             col_center_i++)
        {
            center_signal_type = s4strain[BOARD_STRAIN2(row_center_i,
                                                        col_center_i)];
            for (row_rad_i = row_center_i - S_rad;
                 row_rad_i < row_center_i + S_rad + 1;
                 row_rad_i++)
            {
                for (col_rad_i = col_center_i - S_rad;
                     col_rad_i < col_center_i + S_rad + 1;
                     col_rad_i++)
                {
                    row_rad_i_t = pyModulus(row_rad_i, Nboard_strain[0]);
                    col_rad_i_t = pyModulus(col_rad_i, Nboard_strain[1]);
                    printf("%d\n", BOARD_SIGNAL_NUM3(center_signal_type,
                                                     row_rad_i_t,
                                                     col_rad_i_t) + 1);
                    BOARD_SIGNAL_NUM3(center_signal_type,
                                      row_rad_i_t,
                                      col_rad_i_t) = BOARD_SIGNAL_NUM3(center_signal_type,
                                                                       row_rad_i_t,
                                                                       col_rad_i_t) + 1;
                    printf("%d\n", BOARD_SIGNAL_NUM3(center_signal_type,
                                                     row_rad_i_t,
                                                     col_rad_i_t) + 1);
                }
            }
        }
    }
    '''
    sp.weave.inline(code, ['board_strain', 'board_signal_num', 'S_rad'], support_code=support_code)
                    

    #for row_center_i in range(board_strain.shape[0]):
    #    for col_center_i in range(board_strain.shape[1]):
    #        center_signal_type = s4strain[board_strain[row_center_i,
    #                                                   col_center_i]]
    #        for row_rad_i in range(row_center_i - S_rad, row_center_i + S_rad + 1):
    #            for col_rad_i in range(col_center_i - S_rad, col_center_i + S_rad + 1):
    #                row_rad_i_t = row_rad_i % board_strain.shape[0]
    #                col_rad_i_t = col_rad_i % board_strain.shape[1]
    #                board_signal_num[center_signal_type,
    #                                 row_rad_i_t,
    #                                 col_rad_i_t] = board_signal_num[center_signal_type,
    #                                                                 row_rad_i_t,
    #                                                                 col_rad_i_t] + 1
    return board_signal_num


def init_board_prod(board_strain,
                    board_signal_num,
                    S_th):
    for row_center_i in range(board_strain.shape[0]):
        for col_center_i in range(board_strain.shape[1]):
            receptor_strain = r4strain[board_strain[row_center_i,
                                                    col_center_i]]
            print board_signal_num[receptor_strain,
                                   row_center_i,
                                   col_center_i]
            board_prod[row_center_i,
                       col_center_i] = S_th[board_signal_num[receptor_strain,
                                                             row_center_i,
                                                             col_center_i]]


def init_board_pg_num(board_strain,
                      board_prod,
                      G_th,
                      G_rad):
    for row_center_i in range(board_strain.shape[0]):
        for col_center_i in range(board_strain.shape[1]):
            for row_rad_i in range(row_center_i - G_rad,
                                   row_center_i + G_rad + 1):
                for col_rad_i in range(col_center_i - G_rad,
                                       col_center_i + G_rad + 1):
                    row_rad_i_t = row_rad_i % board_strain.shape[0]
                    col_rad_i_t = col_rad_i % board_strain.shape[1]
                    board_pg_num[row_center_i,
                                 col_center_i] += G_th[board_prod[row_rad_i_t,
                                                                  col_rad_i_t]]


@profile
def init_boards(board_strain,
                board_signal_num,
                board_pg_num,
                board_prod,
                S_th,
                G_th,
                S_rad,
                G_rad):
    init_board_signal_num(board_strain,
                          board_signal_num,
                          S_rad)

    init_board_prod(board_strain,
                    board_signal_num,
                    S_th)

    init_board_pg_num(board_strain,
                      board_prod,
                      G_th,
                      G_rad)

    return (board_signal_num,
            board_prod,
            board_pg_num)


#@numba.jit('i8[:](u8[:],i8,i8)')
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
        'basal_cost':         100, # Basal metabolic cost

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
              'basal_cost': basal_cost,
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
              'board_prod': board_prod,
              'board_pg_num': board_pg_num,
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



# Either load an existing simulation, or create a new one.
if os.path.exists(config['data_filename']):
    # load an existing simulation
    print 'os.path.exists(' + config['data_filename'] + ')'
    data = {key: val for key, val in np.load(config['data_filename']).items()}

    #mt_0 = data['mt_0']
    #mt = data['mt']
    S_cost = data['S_cost']
    R_cost = data['R_cost']
    C_cost = data['C_cost']
    basal_cost = data['basal_cost']
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
    board_strain = data['board_strain']

    board_signal_num = np.zeros((2,
                                 board_size,
                                 board_size),
                                dtype=board_strain.dtype)
    board_pg_num = np.zeros((board_size,
                             board_size),
                            dtype=board_strain.dtype)
    board_prod = np.zeros((board_size,
                           board_size),
                          dtype=board_strain.dtype)
    init_boards(board_strain,
                S_th,
                G_th,
                S_rad,
                G_rad)
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
    np.random.set_state((data['randomstate_current_0'],
                        data['randomstate_current_1'],
                        data['randomstate_current_2'],
                        data['randomstate_current_3'],
                        data['randomstate_current_4']))
else:
    # Create / initiate a new simulation
    print 'os.path.exists(' + config['data_filename'] + ')'
    #########
    # model parameters
    #########

    #######
    # Cost of gene expression
    S_cost = int(config['S_cost'])
    R_cost = int(config['R_cost'])
    C_cost = int(config['C_cost'])

    # B for Baseline, basal, "basic metabolic burden"
    basal_cost = int(config['basal_cost'])

    #####
    # benefit from cooperation. "reward factor" in the article.
    benefit = float(config['benefit'])

    ######
    # mutation per generation
    mutation_rate_r = float(config['mutation_rate_r'])
    mutation_rate_s = float(config['mutation_rate_s'])
    mutation_rate_c = float(config['mutation_rate_c'])

    ## neighbours effects' thresholds
    S_th = np.array([0 if i < int(config['S_th']) else 1 for i in range(10)])
    # quorum threshold
    G_th = np.array([0 if i < int(config['G_th']) else 1 for i in range(10)])
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
    randomstate_start = np.random.get_state()

    # A cell can be Signalling and/or Receptive and/or Cooperative
    R = np.random.rand(board_size, board_size) < float(config['initial_receptives_amount'])
    S = np.random.rand(board_size, board_size) < float(config['initial_signallers_amount'])
    R = np.int8(R)
    S = np.int8(S)
    print R
    print S
    board_strain = R + 2 * S
    print "board_strain", type(board_strain), board_strain.dtype
    print board_strain

    board_signal_num = np.zeros((2,
                                 board_size,
                                 board_size),
                                dtype=board_strain.dtype)
    board_pg_num = np.zeros((board_size,
                             board_size),
                            dtype=board_strain.dtype)
    board_prod = np.zeros((board_size,
                           board_size),
                          dtype=board_strain.dtype)
    init_boards(board_strain,
                board_signal_num,
                board_pg_num,
                board_prod,
                S_th,
                G_th,
                S_rad,
                G_rad)

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
#@numba.jit('i8(uint8[:,:,:],i8,i8)')
#def goods_count(board_strain,
#                pos_row,
#                pos_col):
#    """
#    Counts the number of public goods at a certain position
#    """
##     has fitness effect?:
##         surrounded by enough goods producing cells
##     produces goods?:
##         no receptor and cooperator
##         or
##         surrounded by enough signal producing cells and
##             is receptive and
##             is cooperator
#    goods_sum = 0
#    for g_row_i in range(pos_row - G_rad,
#                         pos_row + G_rad + 1):
#        for g_col_i in range(pos_col - G_rad,
#                             pos_col + G_rad + 1):
#            #print (g_row_i, g_col_i), board[g_row_i%board_size, g_col_i%board_size],
#            #print goods_sum, 
#            # Able to produce public goods
#            if int(board[g_row_i % board_size,
#                     g_col_i % board_size,
#                     COOPERATION]) == 1:
#                #print 'C,',
#                # Isn't receptive
#                if int(board[g_row_i % board_size,
#                         g_col_i % board_size,
#                         RECEPTOR]) == 0:
#                    #print 'R,',
#                    goods_sum += 1
#                # Receptive and signal reaches signal threshold
#                elif int(signal_count(g_row_i, g_col_i, 1)) >= S_th:
#                    #print 'S_th+',
#                    goods_sum += 1
#                #else:
#                    #print 'S_th-',
#            #else:
#                #print 'c',
#            #print goods_sum
#    return goods_sum



#@numba.autojit
#@numba.jit('f8(i8[:,:],i8[:,:],i8[:,:],i8,i8)')
def fitness(board_prod,
            board_pg_num,
            pos_row,
            pos_col):
    pos_row_t = pos_row % board_size
    pos_col_t = pos_col % board_size

    #print 'fitness(),', 'pos_row_t, ', pos_row_t
    #print 'fitness(),', 'pos_col_t, ', pos_col_t

    goods_num = board_prod[pos_row_t,
                           pos_col_t]
    cost = basal_cost
    cost = cost + C_cost * board_prod[pos_row_t,
                                      pos_col_t]
    cost = cost * (1 -  benefit * G_th[goods_num])

    return basal_cost / cost

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

    #result = benefit * (goods_num >= G_th)
    #result = 1 - result
    #result = (basal_cost + result) / result

    ##print 'fitness(),', 'result,', result#, stuffprint(result)
    #return result

def stuffprint(a):
    print 'stuffprint(),', type(a), a.dtype, a.shape


#@numba.jit('i8[:,:](uint8[:,:,:], i8,i8,i8,i8, f8,f8)')
#@numba.jit('i8(uint8[:,:,:],i8,i8,i8,i8,f8,f8)')
def a_competition(board,
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


#@numba.jit('i8(int8[:,:],f8[:,:,:],i8[:])')
def mutate(strain, p_muts, pair_i):
    """
    mutate(board, self, pos) -> NoneType

    For each value of self.board at position "pos",
        change its value at probability of self.mutation_rate_[r/s/c].
    """
    r, s = r4strain[strain], s4strain[strain]
    if int(2 * mutation_rate_r - p_muts[pair_i, 0]):
        r = r4strain[strain] - 1
    if int(2 * mutation_rate_s - p_muts[pair_i, 1]):
        s = r4strain[strain] - 1

    return strain_spec[r,s]

    #pos_row_t = pos_row % board_size
    #pos_col_t = pos_col % board_size
    #strain = board_strain[pos_row, pos_col]
    #r = r4strain[strain]
    #if float(p_r) < float(mutation_rate_r):
    #    if board[pos_row_t, pos_col_t, 0]:
    #        board[pos_row_t, pos_col_t, 0] = 0
    #    else:
    #        board[pos_row_t, pos_col_t, 0] = 1
    #if float(p_s) < float(mutation_rate_s):
    #    if board[pos_row_t, pos_col_t, 1]:
    #        board[pos_row_t, pos_col_t, 1] = 0
    #    else:
    #        board[pos_row_t, pos_col_t, 1] = 1
    #if float(p_c) < float(mutation_rate_c):
    #    if board[pos_row_t, pos_col_t, 2]:
    #        board[pos_row_t, pos_col_t, 2] = 0
    #    else:
    #        board[pos_row_t, pos_col_t, 2] = 1
    ##return board


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
#@numba.jit('i4[:,:](i4[:,:],i4[:])')
def rel_pos_find(a, b, directions):
    """
    Takes two (n, 2) integer arrays (*a* and *b*), and according to a (n,) integers array holding the numbers [0,4] denoting directions.
    """
    code = '''
    int pair_i = 0;
    for (pair_i = 0; pair_i < Na[0]; pair_i++) {
        if (DIRECTIONS1(pair_i) == 0) {
            B2(pair_i, 0) = pyModulus(A2(pair_i, 0) - 1, board_size);
            B2(pair_i, 1) = A2(pair_i, 1);
            continue;
        }
        if (DIRECTIONS1(pair_i) == 1) {
            B2(pair_i, 0) = A2(pair_i, 0);
            B2(pair_i, 1) = pyModulus(A2(pair_i, 1) - 1, board_size);
            continue;
        }
        if (DIRECTIONS1(pair_i) == 2) {
            B2(pair_i, 0) = pyModulus(A2(pair_i, 0) + 1, board_size);
            B2(pair_i, 1) = A2(pair_i, 1);
            continue;
        }
        B2(pair_i, 0) = A2(pair_i, 0);
        B2(pair_i, 1) = pyModulus(A2(pair_i, 1) + 1, board_size);
    }
    '''

    support_code = r'''
    int pyModulus(int a, int b) {
        return ((a % b) + b) % b;
    }
    '''

    sp.weave.inline(code, ['a', 'b', 'directions', 'board_size'], support_code=support_code)

#@numba.autojit
def same(board, pair_i, pos_1s, pos_2s):
    for gene_i in range(3):
        if (board[pos_1s[pair_i, 0], pos_1s[pair_i, 1], gene_i] !=
            board[pos_2s[pair_i, 0], pos_2s[pair_i, 1], gene_i]):
            return 0
    return 1


#@numba.autojit
def update_boards(board_strain,
                  board_signal_num,
                  board_prod,
                  board_pg_num,
                  new_strain,
                  old_strain,
                  pos_row,
                  pos_col):
    """
    Returns the updated signal number, production activity and
    number of public goods boards.
    """
    s4strain = np.array([0, 0, 1, 1], dtype='int8')  # 0: S0, 1: S0, 2: S1, 3: S1 
    for row_center_i in range(pos_row - S_rad, pos_row + S_rad + 1):
        for col_center_i in range(pos_col - S_rad, pos_col + S_rad + 1):
            # t for toroid
            row_center_i_t = row_center_i % board_strain.shape[0]
            col_center_i_t = col_center_i % board_strain.shape[1]

            # Update signal number
            board_signal_num[s4strain[old_strain],
                             row_center_i_t,
                             col_center_i_t] = board_signal_num[s4strain[old_strain],
                                                              row_center_i_t,
                                                              col_center_i_t] - 1
            board_signal_num[s4strain[new_strain],
                             row_center_i_t,
                             col_center_i_t] = board_signal_num[s4strain[new_strain],
                                                              row_center_i_t,
                                                              col_center_i_t] + 1

            ## Update goods production activity
            #old_prod = board_prod[row_center_i_t,
            #                      col_center_i_t]
            #receptor_type = r4strain[board_strain[row_center_i_t,
            #                         col_center_i_t]]
            #new_prod = G_th[board_signal_num[receptor_type,
            #                                 row_center_i_t,
            #                                 col_center_i_t]]
            #board_prod[row_center_i_t,
            #           col_center_i_t] = new_prod

            #if old_prod != new_prod:
            #    pg_difference = new_prod - old_prod
            #    for row_rad_i in range(row_center_i - G_rad,
            #                           row_center_i + G_rad + 1):
            #        for col_rad_i in range(col_center_i - G_rad,
            #                               col_center_i + G_rad + 1):
            #            # t for toroid indexes
            #            row_rad_i_t = row_rad_i % board_strain.shape[0]
            #            col_rad_i_t = col_rad_i % board_strain.shape[1]
            #            new_pg_num = board_pg_num[row_rad_i_t,
            #                                      col_rad_i_t] + pg_difference
            #            board_pg_num[row_rad_i_t, col_rad_i_t] = new_pg_num 
    
    return board_signal_num, board_prod, board_pg_num



#@numba.autojit
@profile
def competitionstep(board_strain,
                    board_signal_num,
                    board_prod,
                    board_pg_num):
    """
    Returns the strain, signal amount, public goods production state and
    public goods amount boards after one generation.

    Runs the competitions and mutation code.
    """
    pos_1s = np.random.randint(board_size, size=(board_size ** 2, 2))
    rel_pos_s = np.random.randint(4, size=(board_size ** 2))
    p_pairs = np.random.rand(board_size ** 2, 2)
    p_muts = np.random.rand(board_size ** 2, 3)
    pos_2s = np.empty(pos_1s.shape, dtype=pos_1s.dtype)
    rel_pos_find(pos_1s, pos_2s, rel_pos_s)
    
    
    for pair_i in range(board_size ** 2):
        #print 'competitionstep(), pair_i first,', pair_i
        a_r = pos_1s[pair_i, 0]
        a_c = pos_1s[pair_i, 1]
        b_r = pos_2s[pair_i, 0]
        b_c = pos_2s[pair_i, 1]
        b_rt= pos_2s[pair_i, 0] % board_strain.shape[0]
        b_ct= pos_2s[pair_i, 1] % board_strain.shape[1]

        #print a_r, a_c
        #print b_r, b_c
        #print board_strain[a_r, a_c]
        #print board_strain[b_rt, b_ct]
        #print same_strain[board_strain[a_r, a_c],
        #                  board_strain[b_rt, b_ct]]

        if same_strain[board_strain[a_r, a_c],
                       board_strain[b_rt, b_ct]]:
            newstrain = mutate(board_strain[a_r,
                                            a_c],
                               p_muts,
                               pair_i)
            board_strain[a_r,
                         a_c] = newstrain
            #for n in [board_strain, board_signal_num, board_prod, board_pg_num, newstrain, board_strain[b_rt, b_ct], a_r, a_c]:
            #    print type(n)
            #    print newstrain, board_strain[b_rt, b_ct], a_r, a_c
            (board_signal_num,
             board_prod,
             board_pg_num) = update_boards(board_strain,
                                           board_signal_num,
                                           board_prod,
                                           board_pg_num,
                                           newstrain,
                                           board_strain[b_rt,
                                                        b_ct],
                                           a_r,
                                           a_c)
            continue

        a_fitness = fitness(board_prod,
                            board_pg_num,
                            a_r,
                            a_c)

        b_fitness = fitness(board_prod,
                            board_pg_num,
                            b_r,
                            b_c)

        a_score = a_fitness * p_pairs[pair_i, 0]
        b_score = b_fitness * p_pairs[pair_i, 1]

        if (float(a_score) > float(b_score)):
            newstrain = mutate(board_strain[b_rt, b_ct],
                               p_muts,
                               pair_i)
            board_strain[b_rt, b_ct] = newstrain
        else:
            newstrain = mutate(board_strain[a_r, a_c],
                               p_muts,
                               pair_i)
            board_strain[a_r, a_c] = newstrain

        return board_strain, board_signal_num, board_prod, board_pg_num

        #equal = sp.weave.inline(r'''
        #                   int gene_i;
        #                   return_val = 1;
        #                   for (gene_i = 0; gene_i < 3; gene_i++) {
        #                       if (BOARD3(a_r, a_c, gene_i) != BOARD3(b_r, b_c, gene_i)) {
        #                           return_val = 0;
        #                           break;
        #                       }
        #                   }''',
        #          ['a_r', 'a_c', 'b_r', 'b_c', 'board'],
        #          {'a_r': int(pos_1s[pair_i, 0]),
        #           'a_c': int(pos_1s[pair_i, 1]),
        #           'b_r': int(pos_2st[pair_i, 0]),
        #           'b_c': int(pos_2st[pair_i, 1]),
        #           'board': board})
        #print 'competitionstep(), pair_i equal,', pair_i
        #print 'competitionstep(),', equal, board[pos_1s[pair_i], pos_1s[pair_i], :], board[pos_2s[pair_i], pos_2s[pair_i], :]
        #if equal:
        #    #print 'competitionstep(), pair_i if equal, ', pair_i
        #    mutate(board,
        #           pos_1s[pair_i, 0],
        #           pos_1s[pair_i, 1],
        #           p_muts[pair_i, 0],
        #           p_muts[pair_i, 1],
        #           p_muts[pair_i, 2])
        #    #print 'competitionstep(), pair_i if equal end, ', pair_i
        #else:
        #    #print 'competitionstep(), pair_i if not equal, ', pair_i
        #    competition(board,
        #                                                                 pos_1s[pair_i, 0],
        #                                                                 pos_1s[pair_i, 1],
        #                                                                 pos_2s[pair_i, 0],
        #                                                                 pos_2s[pair_i, 1],
        #                                                                 p_pairs[pair_i, 0],
        #                                                                 p_pairs[pair_i, 1])
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

def diffuse():
    pass

#@numba.autojit#('(i4[:,:],i4[:,:,:],i4[:,:],i4[:,:],i4[:],i4[:,:])')
@profile
def rotquad90(board,
              directions,
              positions):
    """
    rotquad90(self, direction, position) -> NoneType

    Turns a 2 by 2 sub-array of self.board by 90 degrees.
    Direction:
        0 - turn anticlockwise
        1 - turn clockwise
    with its lowest valued coordinate (upper left) in the "position" coordinate value.
    """
    #l = [board_strain,
    #          board_signal_num,
    #          board_prod,
    #          board_pg_num,
    #          directions,
    #          positions]
    #for i in range(len(l)):
    #    print type(l[i]), l[i].shape, l[i].dtype
    code = r'''
KEEP_VALS2(0,0) = BOARD_STRAIN2((int) row0,(int)  col0);
KEEP_VALS2(0,1) = BOARD_STRAIN2((int) row0,(int)  col1);
KEEP_VALS2(1,0) = BOARD_STRAIN2((int) row1,(int)  col0);
KEEP_VALS2(1,1) = BOARD_STRAIN2((int) row1,(int)  col1);
		
if ((int) clockwise) {
    // clockwise
    // # 00 01 -> 01 11
    // # 10 11 -> 00 10
    NEW_VALS2(0,0) = KEEP_VALS2(0,1);
    NEW_VALS2(0,1) = KEEP_VALS2(1,1);
    NEW_VALS2(1,1) = KEEP_VALS2(1,0);
    NEW_VALS2(1,0) = KEEP_VALS2(0,0);
} else {
    // anticlockwise
    // # 00 01 -> 10 00
    // # 10 11 -> 11 01
    NEW_VALS2(0,0) = KEEP_VALS2(1,0);
    NEW_VALS2(1,0) = KEEP_VALS2(1,1);		
    NEW_VALS2(1,1) = KEEP_VALS2(0,1);
    NEW_VALS2(0,1) = KEEP_VALS2(0,0);		
}
		
BOARD_STRAIN2((int) row0, (int) col0) = NEW_VALS2(0,0);
BOARD_STRAIN2((int) row0, (int) col1) = NEW_VALS2(0,1);
BOARD_STRAIN2((int) row1, (int) col0) = NEW_VALS2(1,0);
BOARD_STRAIN2((int) row1, (int) col1) = NEW_VALS2(1,1);
/*
        update_boards(new py::tuple(py_board_strain,
                      py_board_signal_num,
                      py_board_prod,
                      py_board_pg_num,
                      NEW_VALS2(0,0),
                      KEEP_VALS2(0,0),
                      row0,
                      col0))

        update_boards(py_board_strain,
                      py_board_signal_num,
                      py_board_prod,
                      py_board_pg_num,
                      NEW_VALS2(0,1),
                      KEEP_VALS2(0,1),
                      row0,
                      col1)

        update_boards(py_board_strain,
                      py_board_signal_num,
                      py_board_prod,
                      py_board_pg_num,
                      NEW_VALS2(1,0),
                      KEEP_VALS2(1,0),
                      row1,
                      col0)

        update_boards(py_board_strain,
                      py_board_signal_num,
                      py_board_prod,
                      py_board_pg_num,
                      NEW_VALS2(1,1),
                      KEEP_VALS2(1,1),
                      row1,
                      col1)
                      */
    '''
    new_vals = np.empty((2,2), dtype=board_strain.dtype)
    keep_vals = np.empty((2,2), dtype=board_strain.dtype)
    for rotation_i in range(diffusion_step_num):
        row0 = sp.weave.inline('return_val = POSITIONS2(rotation_i, axis);', ['positions', 'rotation_i','axis'],{'rotation_i':rotation_i,
                                                                                                       'positions':positions,
                                                                                                       'axis': 0})
        col0 = sp.weave.inline('return_val = POSITIONS2(rotation_i, axis);', ['positions', 'rotation_i','axis'],{'rotation_i':rotation_i,
                                                                                                       'positions':positions,
                                                                                                       'axis': 0})
        row1 = sp.weave.inline('return_val = ((a + 1) % b + b) % b;', ['a','b'],{'a':row0,'b': board_strain.shape[0]})
        col1 = sp.weave.inline('return_val = ((a + 1) % b + b) % b;', ['a','b'],{'a':row0,'b': board_strain.shape[0]})
        sp.weave.inline(code,
                        ['board_strain',
                         'board_signal_num',
                         'board_prod',
                         'board_pg_num',
                         'clockwise',
                         'new_vals',
                         'keep_vals',
                         'row0',
                         'row1',
                         'col0',
                         'col1',
                         'diffusion_step_num'],
                        {'row0': row0,
                         'row1': row1,
                         'col0': col0,
                         'col1': col1,
                         'board_strain': board_strain,
                         'board_signal_num': board_signal_num,
                         'board_prod': board_prod,
                         'board_pg_num': board_pg_num,
                         'clockwise': directions[rotation_i],
                         'new_vals': new_vals,
                         'keep_vals': keep_vals,
                         'diffusion_step_num': diffusion_step_num})
    #    (board_signal_num, board_prod, board_pg_num) = update_boards(board_strain, board_signal_num, board_prod, board_pg_num, new_vals[0,0], keep_vals[0,0], row0, col0)

    #    (board_signal_num,
    #     board_prod,
    #     board_pg_num) = update_boards(board_strain,
    #                  board_signal_num,
    #                  board_prod,
    #                  board_pg_num,
    #                  new_vals[0,1],
    #                  keep_vals[0,1],
    #                  row0,
    #                  col1)

    #    (board_signal_num,
    #     board_prod,
    #     board_pg_num) = update_boards(board_strain,
    #                  board_signal_num,
    #                  board_prod,
    #                  board_pg_num,
    #                  new_vals[1,0],
    #                  keep_vals[1,0],
    #                  row1,
    #                  col0)

    #    (board_signal_num,
    #     board_prod,
    #     board_pg_num) = update_boards(board_strain,
    #                  board_signal_num,
    #                  board_prod,
    #                  board_pg_num,
    #                  new_vals[1,1],
    #                  keep_vals[1,1],
    #                  row1,
    #                  col1)
    #l = [board_strain,
    #          board_signal_num,
    #          board_prod,
    #          board_pg_num,
    #          directions,
    #          positions]
    #for i in range(len(l)):
    #    print type(l[i]), l[i].shape, l[i].dtype

    return (board_strain, board_signal_num, board_prod, board_pg_num)

    #if direction < 0.5:
    #    board[row0, col0, :], board[row0, col1, :], board[row1, col0, :], board[row1, col1, :] = board[row0, col1, :], board[row1, col1, :], board[row0, col0, :], board[row1, col0, :]
    #                                          # 00 01 -> 01 11
    #                                          # 10 11 -> 00 10
    #else:
    #    board[row0, col0, :], board[row0, col1, :], board[row1, col0, :], board[row1, col1, :] = board[row1, col0, :], board[row0, col0, :], board[row1, col1, :], board[row0, col1, :]
    #                                          # 00 01 -> 10 00
    #                                          # 10 11 -> 11 01


#@numba.autojit#('void(i4[:,:])') # , uint64[:])')
@profile
def nextstep(board_strain,
             board_signal_num,
             board_prod,
             board_pg_num):
    competition_result = competitionstep(board_strain,
                                     board_signal_num,
                                     board_prod,
                                     board_pg_num)
    print type(competition_result
              )
    board_strain = competition_result[0]
    board_signal_num = competition_result[1]
    board_prod = competition_result[2]
    board_pg_num = competition_result[3]

    directions = np.random.randint(2, size=diffusion_step_num)
    positions = np.random.randint(board_size, size=(diffusion_step_num, 2))

    (board_strain,
     board_signal_num,
     board_prod,
     board_pg_num) = rotquad90(board_strain,
                               board_signal_num,
                               board_prod,
                               board_pg_num,
                               directions,
                               positions)

    return (board_strain,
            board_signal_num,
            board_prod,
            board_pg_num)

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


#@numba.autojit#('void(int8[:,:,:],f8[:])')
@profile
def go(board_strain,
       board_signal_num,
       board_prod,
       board_pg_num,
       S_th,
       S_rad,
       G_th,
       G_rad,
       times):
    print 'start go()'
    #    every = 30*60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    print 'go(),', board_strain.shape[0], generations
    times[0] = np.float64(time.time())
    for step_count in range(generations):
        #(board_strain,
        # board_signal_num,
        # board_prod,
        # board_pg_num) = 
        nextstep(board_strain,
                                  board_signal_num,
                                  board_prod,
                                  board_pg_num)
        #print board_strain
        times[step_count + 1] = time.time()
        if (step_count % 5) == 0:
            print 'go(),step:', step_count
            print 'go(),300', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 300**2 / 60 / 60
            print 'go(),300', np.average(times[2:step_count + 2] - times[1:step_count+1]) * 10000 / board_size**2 * 300**2 / 60 / 60
            print 'go(),150', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 150**2 / 60 / 60
            print 'go(),150', np.average(times[2:step_count + 2] - times[1:step_count+1]) * 10000 / board_size**2 * 150**2 / 60 / 60
            print 'go(),75', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 75**2 / 60 / 60
            print 'go(),75', np.average(times[2:step_count + 2] - times[1:step_count+1]) * 10000 / board_size**2 * 75**2 / 60 / 60
            print 'go(),50', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 50**2 / 60 / 60
            print 'go(),50', np.average(times[2:step_count + 2] - times[1:step_count+1]) * 10000 / board_size**2 * 50**2 / 60 / 60
            print 'go(),25', (times[step_count + 1] - times[step_count]) * 10000 / board_size**2 * 25**2 / 60 / 60
            print 'go(),25', np.average(times[2:step_count + 2] - times[1:step_count+1]) * 10000 / board_size**2 * 25**2 / 60 / 60
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
go(board_strain,
   board_signal_num,
   board_prod,
   board_pg_num,
   S_th,
   S_rad,
   G_th,
   G_rad,
   times)
def timetofinish(timepoints, board_size, generations):
    np.average(timepoints[1:] - timepoints[:-1])

l=[]
for time_0, time_1 in zip(times[:-1], times[1:]):
    l.append(time_1-time_0)
print 'root level,', np.average(l[2:]) * 10000 / board_strain.shape[0]**2 * 300**2 / 60 / 60 

    
np.savez(config['data_filename'], **get_state())
