#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import h5py
#import os
import time
import numpy as np
import scipy as sp
import scipy.signal
import scipy.weave
import pylab as pl
#import timeit
import sys
import signal
import ConfigParser

labels = ['Ignorant (csr)', 'Voyeur (csR)', 'Liar (cSr)', 'Lame (cSR)',
          'Blunt (Csr)', 'Shy (CsR)', 'Vain (CSr)', 'Honest (CSR)']


class Strife:

    """
    Strife class

    When initiated with no config dictionary, default_config is loaded.
    """

    def __init__(self, config=None):
        """
        Creates a Game of Strife board. A dictionary argument called "config".
        """
        if config is None:
            config = default_config

        ########
        # filname
        self.data_filename = sp.array(config['data_filename'])

        #########
        # model parameters
        #########

        #######
        # Cost of gene expression
        self.S_cost = sp.float64(config['S_cost'])
        self.R_cost = sp.float64(config['R_cost'])
        self.C_cost = sp.float64(config['C_cost'])
        self.metabolic_baseline = sp.float64(config['metabolic_baseline'])  # B for Baseline, basal, "basic metabolic burden"

        #####
        # benefit from cooperation. "reward factor" in the article.
        self.benefit = sp.float64(config['benefit'])

        ######
        # mutation per generation
        self.mutation_rate_r = sp.float64(config['mutation_rate_r'])
        self.mutation_rate_s = sp.float64(config['mutation_rate_s'])
        self.mutation_rate_c = sp.float64(config['mutation_rate_c'])

        ## neighbours effects' thresholds
        self.S_th = sp.int64(config['S_th'])
        # quorum threshold
        self.C_th = sp.int64(config['C_th'])
        # Cooperation threshold. Above it, public goods makes a difference.

        ## Probability of each single diffusion operation
        self.diffusion_amount = sp.float64(config['diffusion_amount'])

        #######
        # radius of Signal or Cooperation effects.
        self.S_rad = sp.int64(config['S_rad'])
        self.C_rad = sp.int64(config['C_rad'])

#        #######
#        # We don't want to create these again and again while counting for each competition.
#        self.horizontal_s_count = sp.empty(2 + 2 * self.S_rad + 2 * self.C_rad,
#                                           1 + 2 * self.S_rad + 2 * self.C_rad)
#        self.vertical_s_count = sp.empty(2 + 2 * self.S_rad + 2 * self.C_rad,
#                                         1 + 2 * self.S_rad + 2 * self.C_rad)
#        self.horizontal_c_count = sp.empty(2 + 2 * self.C_rad,
#                                           1 + 2 * self.C_rad)
#        self.vertical_c_count = sp.empty(2 + 2 * self.C_rad,
#                                         1 + 2 * self.C_rad)

        self.generations = sp.int64(config['generations'])

        self.board_size = sp.int64(config['board_size'])

        #####
        # settings
        #####

        self.NEIGHBOUR_REL_POS = sp.array([(0, -1), (0, 1), (-1, 0), (1, 0)])

        ## time keeping
        # number of generations the simulation will run
        # each generation is defined as the average number of steps for which
        # each cell on the board was in a competition once since last
        # generation (?)

        #######
        # We'll increase step_count by one at the end of each step function.
        self.step_count = sp.int64(0)

        self.steps_final = sp.int64(self.generations * self.board_size ** 2)

        # A cell can be Signalling and/or Receptive and/or Cooperative
        R = sp.rand(self.board_size, self.board_size) < config['initial_receptives_amount']
        S = sp.rand(self.board_size, self.board_size) < config['initial_signallers_amount']
        C = sp.rand(self.board_size, self.board_size) < config['initial_cooperators_amount']
        self.board = sp.array([R, S, C]).transpose((1, 2, 0))
        assert self.board_size == self.board.shape[0] and \
               self.board_size == self.board.shape[1], \
            'B.shape: {0}\nN: {1}\nWanted: {2}'.format(self.board.shape,
                                                       self.board_size,
                                                       (self.board_size,
                                                        self.board_size,
                                                        3))

        self.genotype_num = sp.int64(8)

        ## data sampling
        # we will take a frequency sample some number of times per generation
        self.steps_per_gen = sp.int64(self.board_size ** 2)
        self.samples_per_gen = sp.int64(1)
        self.samples_num = sp.int64(self.samples_per_gen * self.generations)
        self.samples_board_num = sp.int64(self.samples_num // 10)
        self.steps_per_sample = sp.int64(sp.floor(1.0 * self.steps_per_gen // self.samples_per_gen), dtype=sp.int64)
        self.steps_per_board_sample = sp.int64(10 * self.steps_per_sample)
        self.sample_count = sp.int64(0)
        # We want to know the frequency of each genotype per generation
        self.samples_frequency = sp.empty((self.samples_num, self.genotype_num), dtype='int64')
        self.samples_nhood = sp.empty((self.samples_num, self.genotype_num, self.genotype_num), dtype=sp.int64)
        self.samples_board = sp.empty((self.samples_board_num, self.board_size, self.board_size, 3), dtype=sp.int64)

        support_code = open('support_code.c').read().format({'S_rad': self.S_rad,
                                                             'C_rad': self.C_rad,
                                                             'Nboard_row': board.shape[0],
                                                             'Nboard_col': board.shape[1],
                                                             'RECEPTOR': 0,
                                                             'SIGNAL': 1,
                                                             'COOPERATION': 2,
                                                             'R_cost': self.R_cost,
                                                             'S_cost': self.S_cost,
                                                             'C_cost': self.C_cost})

        code = '''
long C_th = {C_th};
long benefit = {benefit};
long metabolism_1 = (c_count(C_POS1(0), C_POS1(1)) >= C_th) * (1 - (benefit);
long metabolism_2 = ();
if ((metabolism(C_POS_11[0], C_POS_11[1]) * p1) >
    (metabolism(C_POS_11[0], C_POS_11[1]) * p1))
    // The number at the end of the "C_POS" denotes number dimensions in the nd-array (ndim).
{
    return 1;
}
else
{
    return 2;
}'''

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
        assert p1.shape == () and p2.shape == (), \
            'p1 ({p1}) and p2 ({p2}) need to be of shape ().\
Are actually of {p1shape} and {p2shape}'.format(p1=p1,
                                                p2=p2,
                                                p1shape=p1.shape,
                                                p2shape=p2.shape)
        assert p1.dtype.kind == 'f' and p2.dtype.kind == 'f', \
            'p1 ({0}) and p2 ({1}) need to be of float dtype,\
but are actually of {p1dtype} and {p2dtype}.'.format(p1, p2, p1.dtype, p2.dtype)
        assert ((0 <= p1) & (p1 < 1) and (0 <= p2) & (p2 < 1)).all(), 'p1 ({0}) and p2 ({1}) need to be over [0, 1)'.format(p1, p2)

        # c_pos_2's coordinates in a torus:
        c_pos_2t = c_pos_2 % self.board_size
        assert (0 <= c_pos_1[0]) and \
               (0 <= c_pos_1[1]) and \
               (0 <= c_pos_2t[0]) and \
               (0 <= c_pos_2t[1]) and \
               (c_pos_1[0] < self.board_size) and \
               (c_pos_1[1] < self.board_size) and \
               (c_pos_2t[0] < self.board_size) and \
               (c_pos_2t[1] < self.board_size), 'c_pos_1: {0}\nc_pos_2t: {1}'.format(c_pos_1, c_pos_2t)

        winner = scipy.weave.inline(code, ['board',
                                           'row',
                                           'col',
                                           'c_pos_1',
                                           'c_pos_2',
                                           'B_cost',
                                           'R_cost',
                                           'S_cost',
                                           'C_cost'],
                                    {'board': self.board,
                                     'row': row,
                                     'col': col,
                                     'c_pos_1': c_pos_1,
                                     'c_pos_2': c_pos_2,
                                     'B_cost': self.B_cost,
                                     'R_cost': self.R_cost,
                                     'S_cost': self.S_cost,
                                     'C_cost': self.C_cost},
                                    support_code=support_code)

        return (c_pos_1, c_pos_2) if (winner == 1) else
               (c_pos_2, c_pos_1)

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

    def save_h5(self):
        """
        Saves the attributes of self whose names show up as keys in self.parameters.
        """
        with h5py.File(self.data_filename) as ff:
            for key, val in vars(self).items():
                if type(val).__module__ == np.__name__:
                    try:
                        print(key, type(getattr(self, key)))
                        ff[key] = getattr(self, key)
                    except:
                        ff[key][...] = getattr(self, key)

    def load_h5(self):
        with h5py.File(self.data_filename) as ff:
            d = {key: val[...] for key, val in ff.items()}
            self.__init_unpack_parameters__(d)

    def sample(self):
        joint_board = self.board[:, :, 0] + 2 * self.board[:, :, 1] + 4 * self.board[:, :, 2]
        for gene in range(8):
            gene_board = joint_board == gene
            gene_frequency = sp.sum(gene_board)
            # neighbours_genotype
            for nh_gene in range(8):
                nh_board = joint_board == nh_gene
                nh_gene_count = sp.signal.convolve2d(nh_board, sp.ones((3, 3)),
                                                     mode='same', boundary='wrap')
                nh_gene_count_of_gene = sp.sum(nh_gene_count * gene_board, dtype=sp.int64)
                self.samples_nhood[self.sample_count, gene, nh_gene] = nh_gene_count_of_gene
            self.samples_frequency[self.sample_count, gene] = gene_frequency
        self.sample_count += 1

    def nextstep(self):
        print('generation:', self.step_count)
        for i in range(self.board_size ** 2):
            print('competition:', i)
            ##
            # Draw two adjacent positions.
            # We'll use relative positions to compute exact positions of 2nd competitor cell
            pos1 = sp.random.randint(self.board_size, size=2)
            pos2 = pos1 + self.NEIGHBOUR_REL_POS[sp.random.randint(4)]
            p_pair = sp.rand(2)
            winner_pos, loser_pos = self.competition(pos1, pos2, p_pair)
            self.copycell(winner_pos, loser_pos)
            self.mutate(loser_pos)
        for i in range(sp.int64((self.board_size ** 2) * self.diffusion_amount // 4)):
            print('diffusion: {0}'.format(i))
            direction = sp.rand()
            position = sp.random.randint(self.board_size, size=2)
            rotquad90(self.board, direction, position)
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
            res[:, i] = sp.array([self.samples_frequency[:, i] + sp.sum(self.samples_frequency[:, :i])])
        return res

    def imagify_data(self):
        ## package boards' data into a displayable array.
        return sp.array([self.S, self.R, self.C])

    def display_frequency_timeseries(self):
        for i in range(8):
            pl.plot(sp.arange(self.samples_num), self.samples_frequency[:, i], label=str(i), fillstyle='bottom')


def go(a):
    signal.signal(signal.SIGINT, handler_maker(a))
    t = time.time()
    every = 30 * 60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
    print("t: {0:f}, steps thus far: {1:d}".format(t, a.step_count))
    steps_a = a.step_count
    while a.step_count <= a.generations:
        a.nextstep()
        delta_t = time.time() - t
        if delta_t > every:
            t = time.time()
            a.save_h5()
            steps_delta = a.step_count - steps_a
            steps_a = a.step_count
            eta = 1.0 * delta_t / (steps_delta + 1) * (a.steps_final - a.step_count)
            print("t: {0:f}, approx. time to fin: {1:f}".format(t, eta))
            print("steps taken = {0}, steps since last save = {1}".format(a.step_count, steps_delta))
            sys.exit(1)


# TODO: Handler of signals.
def handler_maker(a_game):
    def handler(signum, frame):
        print('Signal handler called with signal', signum)
        a_game.save_h5()
        print('game saved')
        raise
    return handler


def assert_ndim(arr, nd):
    assert arr.ndim == nd, 'Wrong number of dimensions.\nExists {0} but {1} is wanted.'.format(arr.ndim, nd)


def assert_shape(arr, shape):
    assert arr.shape == shape, 'Wrong shape\nExists {0} but {1} is wanted.'.format(arr.shape, shape)

default_config = {
        'S_cost':                1,    # Metabolic cost of signalling
        'R_cost':                3,    # Metabolic cost of having a receptor
        'C_cost':                30,   # Metabolic cost of being cooperative
        'metabolic_baseline':    100,  # Basal metabolic cost

        # The fraction reduced, out of total metabolic cost, when
        #    public goods reach threshold of benefit.
        'benefit':               0.9,

        # Likelihoods of switch (on to off and vica versa) for each gene per cell per generation.
        'mutation_rate_r': 1e-4,
        'mutation_rate_s': 1e-4,
        'mutation_rate_c': 1e-4,

        # Amount of signal needed for a receptive and cooperative cell to
        #    start producing public goods.
        'S_th':               3,
        # Amount of public goods needed for metabolic benefit.
        'C_th':               3,

        'diffusion_amount': 0.5,       # Average fraction of cells out of total cells on the board (board_size**2)
                                       #    which will be involved
        'board_size':        10,       # The length of the board. The board will be of shape (board_size, board_size).
        'generations':       10,       # Each generation involves, on average, all cells in a competition, meaning
                                       #    board_size**2 competitions.
        'S_rad':              1,       # Radius of the signal's effect.
        'C_rad':              1,
        'samples_per_gen':    1,
        'initial_receptives_amount': 0,
        'initial_signallers_amount': 0,
        'initial_cooperators_amount': 0,
        'data_filename': 'strife.h5'}


def load_config(config_filename):
    """
    Takes a string holding filename and returns a dictionary with all the
      configuration values.
    """
    our_config = default_config

    config = ConfigParser.SafeConfigParser()
    config.read(config_filename)

    for key, val in our_config:
        our_config[key] = config.get('Config', key)

    return our_config


def rotquad90(board, direction, position):
    """
    rotquad90(self, direction, position) -> NoneType

    Turns a 2 by 2 sub-array of self.board by 90 degrees.
    Direction:
        0 - turn anticlockwise
        1 - turn clockwise
    with its lowest valued coordinate (upper left) in the "position" coordinate
      value.
    """
    inline_code = r'''
long int temp_value;

long int row_i, col_i, gene_i;
row_i = col_i = gene_i = 0;

long int row0, col0;
row0 = POSITION1(0);
col0 = POSITION1(1);

long int row1 = (row0 + 1) % Nboard[0];
long int col1 = (row0 + 1) % Nboard[1];
/*
for (int row_i = 0; row_i < 2; row_i++)
{
    for (int col_i = 0; col_i < 2; col_i++)
    {
        for (int gene_i = 0; gene_i < 3; gene_i++)
        {
        printf("%d", BOARD3((row_i+row0) % Nboard[0],
                            (col_i+col0) % Nboard[1],
                            gene_i));
        }
    printf(" ");
    }
    printf("\n");
}
*/
if (direction == 0)
{
    for (gene_i = 0; gene_i < 3; gene_i++)
    {
        temp_value = BOARD3(row0, col0, gene_i);                   // A  B
                                                                   // C  D

        BOARD3(row0, col0, gene_i) = BOARD3(row0, col1, gene_i);   // B  B
                                                                   // C  D

        BOARD3(row0, col1, gene_i) = BOARD3(row1, col1, gene_i);   // B  D
                                                                   // C  D

        BOARD3(row1, col1, gene_i) = BOARD3(row1, col0, gene_i);   // B  D
                                                                   // C  C

        BOARD3(row1, col0, gene_i) = temp_value;                   // B  C
                                                                   // A  D
    }
}
else
{
    for (gene_i = 0; gene_i < 3; gene_i++)
    {
        temp_value = BOARD3(row0, col0, gene_i);                 // A  B
                                                                 // C  D

        BOARD3(row0, col0, gene_i) = BOARD3(row1, col0, gene_i); // C  B
                                                                 // C  D

        BOARD3(row1, col0, gene_i) = BOARD3(row1, col1, gene_i); // C  B
                                                                 // D  D

        BOARD3(row1, col1, gene_i) = BOARD3(row0, col1, gene_i); // C  B
                                                                 // D  B

        BOARD3(row0, col1, gene_i) = temp_value;                 // C  A
                                                                 // D  B
    }
}
/*
for (int row_i = 0; row_i < 2; row_i++)
{
    for (int col_i = 0; col_i < 2; col_i++)
    {
        for (int gene_i = 0; gene_i < 3; gene_i++)
        {
            printf("%d", BOARD3((row0 + row_i) % Nboard[0],
                                (col0 + col_i) % Nboard[1],
                                gene_i));
        }
        printf(" ");
    }
    printf("\n");
}
*/
'''
    sp.weave.inline(inline_code, ['board',
                                  'position',
                                  'direction',
                                  'row_board_size',
                                  'col_board_size'],
                                 {'board': board,
                                  'position': position,
                                  'direction': direction,
                                  'row_board_size': board.shape[0] ** 3,
                                  'col_board_size': board.shape[1] ** 2})
