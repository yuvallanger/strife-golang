#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
A model by Dr. Avigdor Eldar based on Czárán's work.
http://www.plosone.org/article/info:doi/10.1371/journal.pone.0006655
"""

import h5py
#import os
import time
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
        self.S_cost = sp.int64(config['S_cost'])
        self.R_cost = sp.int64(config['R_cost'])
        self.C_cost = sp.int64(config['C_cost'])
        self.metabolic_baseline = sp.int64(config['metabolic_baseline'])  # B for Baseline, basal, "basic metabolic burden"

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

        #######
        # We don't want to create these again and again while counting for each competition.
        self.horizontal_s_count = sp.empty(2 + 2 * self.S_rad + 2 * self.C_rad,
                                           1 + 2 * self.S_rad + 2 * self.C_rad)
        self.vertical_s_count = sp.empty(2 + 2 * self.S_rad + 2 * self.C_rad,
                                         1 + 2 * self.S_rad + 2 * self.C_rad)
        self.horizontal_c_count = sp.empty(2 + 2 * self.C_rad,
                                           1 + 2 * self.C_rad)
        self.vertical_c_count = sp.empty(2 + 2 * self.C_rad,
                                         1 + 2 * self.C_rad)

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

        # diameter of the convolution matrix
        diameter = lambda x: sp.int64(2 * x + 1)
        S_len = diameter(self.S_rad)
        C_len = diameter(self.C_rad)

        # the convolution matrix used to count neighbours
        self.S_kernel = sp.ones((S_len, S_len))
        self.C_kernel = sp.ones((C_len, C_len))

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
    def count_neighbors(self, rows, cols, gene, allele, radius):
        """
        Counts all neighbors that have "allele" as their "gene" at a Moore's "radius" around each cell that lies between
          rows[0] and rows[1] and between cols[0] and cols[1], inclusive.

        Acts like scipy.signal.convolve2d with mode = 'valid', except that it only counts, does not convolves.
        """
        count_neighbors_code = r'''
long nh_row_i, nh_col_i;
long count = 0;
long nh_row_i_torus; // nh stands for neighborhood
long nh_col_i_torus;
long row_center, col_center;
long rl, rh, cl, ch;


if (C_POS_A1(0) < C_POS_B1(0))
{
    rl = C_POS_A1(0);
    rh = C_POS_B1(0);
}
else
{
    rl = C_POS_B1(0);
    rh = C_POS_A1(0);
}

if (C_POS_A1(0) < C_POS_B1(0))
{
    cl = C_POS_A1(1);
    ch = C_POS_B1(1);
}
else
{
    cl = C_POS_B1(1);
    ch = C_POS_A1(1);
}

/*
        board = array([0, 1, 2, 3, 4])
        I want sums of 3, 4 and 5 so
        The range of cells we want to have sums of neighbors is denoted in "cols".
        Lower is inclusive. Higher is exclusive.
            cols = [3 6)
        The integers we span in the outer loop are exactly the ones in the range of "cols".
            center = [3,6)
        We will span the radius around the centeral cell.
        From "center - radius" inclusive to "center + radius + 1" exclusive.
            col_i = [c-r, c+r+1)
        The board we want to hold the sums must have the shape of the ranges we sample.
        So, the shape in this example must hold 3, 4, and 5.
            cols[1] - cols[0] == 3
        Hence,
            sum_board.shape == (cols[1]-cols[0], )
        We can only get and set our main game board within its range.
            board.shape = (5, )
        This will yield an error.
            BOARD1[-1] = 1
        And so will this.
            BOARD1[5] = 1
        Our game board is torus shaped. In theory, this is always true:
            BOARD1[k * Nboard[0] + n] == BOARD1[n]
        But we live in C, so we must torusify our indexes like so:
            col_i_torus = col_i % Nboard[0] + Nboard[0]
        When we check for equivalence, we use it like this:
            BOARD1[col_i_torus] == foo
        And add one whenever it's true
            count++
        For each col_center, there is one sum_board index that is:
            sum_board_i = col_center - COLS1(0)
        In our example, the three indexes will be:
            3 - 3 == 0
            4 - 3 == 1
            5 - 3 == 2
        So when we're out of the inner loop, we assign to the sum_board.
            SUM_BOARD1(col_center - COLS1(0)) = count

*/
for (row_center = ROWS1(0); row_center < ROWS1(1); row_center++)
{
    for (col_center = COLS1(0); col_center < COLS1(1); col_center++)
    {
        count = 0;
        for (nh_row_i = row_center - radius; nh_row_i < row_center + radius + 1; nh_row_i++)
        {
            for (nh_col_i = col_center - radius; nh_col_i < col_center + radius + 1; nh_col_i++)
            {
                // This is a funny thing. Python's modulus operation always yields a positive number,
                //   both for a positive and a negative first argument.
                //   C/C++'s modulus operation will yield a negative number for a negative first argument.
                //   This is remedied by adding the second argument
                nh_row_i_torus = (nh_row_i < 0) ?
                                 nh_row_i % Nboard[0] + Nboard[0] :
                                 nh_row_i % Nboard[0];
                nh_col_i_torus = (nh_col_i < 0) ?
                                 nh_col_i % Nboard[1] + Nboard[1] :
                                 nh_col_i % Nboard[1];
                printf("nh_row_i: %d, ", nh_row_i);
                printf("nh_row_i_torus: %d, ", nh_row_i_torus);
                printf("gene: %d, ", gene);
                printf("board: %d, ", BOARD3(nh_row_i_torus, nh_col_i_torus, gene));
                printf("allele: %d, ", allele);
                printf("\n");
                if (BOARD3(nh_row_i_torus, nh_col_i_torus, gene) == allele)
                {
                    printf("BOARD3(%d, %d, %d) == %d; allale == %d\n",
                           nh_row_i_torus,
                           nh_col_i_torus,
                           gene,
                           BOARD3(nh_row_i_torus, nh_col_i_torus, gene),
                           allele);
                    count++;
                }
            }
        }
        SUM_BOARD2(row_center - ROWS1(0), col_center - COLS1(0)) = count;
    }
}

// We find public goods producing bacteria.
for (row_center = ROWS1(0); row_center < ROWS1(1); row_center++)
{
    for (col_center = COLS1(0); col_center < COLS1(1); col_center++)
    {
        for (row_i = 0; row_i < 1 + C_rad + row_relation; row_i++)
        {
            for (col_i = 0; col_i < 1 + C_rad + col_relation; col_i)
            {
                // FIXME - The indexing is off
                if ((!BOARD3(row_i, col_i, RECEPTOR) && BOARD3(row_i, col_i, COOPERATION)) ||
                    (BOARD3(row_i, col_i, RECEPTOR) && SUM_BOARD2(row_i, col_i) && BOARD3(row_i, col_i, COOPERATION)))
                {
                    C_BOARD2(row_i, col_i) = 1;
                }
                else
                {
                    C_BOARD2(row_i, col_i) = 0;
                }
            }
        }
    }
}

// We'll count public goods and use it to see
//   if there is a fitness benefit
count = 0;
for (row_center = 0; row_center < 
for (row_i = 0; row_i < 1 + row_relation; row_i++)
{
    for (col_i = 0; col_i < 1 + col_relation; col_i++)
    {
        if (C_BOARD2)
        {
            count++;
        }
    }

    count = 0;
}
        
'''
        sum_board = scipy.empty((rows[1] - rows[0], cols[1] - cols[0]), dtype=sp.int64)
        print(sum_board.shape)
        sp.weave.inline(count_neighbors_code,
                        ['rows', 'cols', 'gene', 'allele', 'board', 'radius', 'sum_board'],
                        {'rows': rows, 'cols': cols, 'gene': gene, 'allele': allele,
                         'board': self.board, 'sum_board': sum_board, 'radius': radius})
        return sum_board

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

        # two identical cells competing will result in two identical cells,
        # so we will return now with no further calculation of this competition.
        if (self.board[c_pos_1[0], c_pos_1[1]] == self.board[c_pos_2t[0], c_pos_2t[1]]).all():
            return (c_pos_2t, c_pos_1)

        ## We will optimize by taking a sub array from each genotype array around the competitors.

        # rl, ch - row low, col high
#        twosort = lambda x, y: (x, y) if x < y else (y, x)
        rl, rh = (c_pos_1[0], c_pos_2[0]) if c_pos_1[0] < c_pos_2[0] else \
                 (c_pos_2[0], c_pos_1[0])
        cl, ch = (c_pos_1[1], c_pos_2[1]) if c_pos_1[1] < c_pos_2[1] else \
                 (c_pos_2[1], c_pos_1[1])

        # For signallers, we take both S_rad and C_rad around our competitors because
        #   signallers affect CR[Ss] cells which, with their public goods, affect our competitors
        s_rows = (rl - self.C_rad - self.S_rad, rh + self.C_rad + self.S_rad + 1)
        s_cols = (cl - self.C_rad - self.S_rad, ch + self.C_rad + self.S_rad + 1)
        c_rows = (rl - self.C_rad, rh + self.C_rad + 1)
        c_cols = (cl - self.C_rad, ch + self.C_rad + 1)

        # I'll use this to keep the arrays as 2d (ndim=2)
        assert rl == rh or cl == ch, 'rl: {0}\nrh: {1}\ncl: {2}\n ch: {3}'.format(rl, rh, cl, ch)

        if rl == rh:      # Competitors are on the same row
            shape = (1, 2)
            s_count = self.horizontal_s_count
            c_count = self.horizontal_c_count
        else:             # otherwise, they're on the same column.
            shape = (2, 1)
            s_count = self.vertical_s_count
            c_count = self.vertical_c_count

        assert self.S_rad.dtype.kind == sp.int8(1).dtype.kind, 'Got {0}, wanted {1}'.format(self.S_rad.dtype.kind, sp.int_(1).dtype.kind)
        assert self.C_rad.dtype.kind == sp.int8(1).dtype.kind, 'Got {0}, wanted {1}'.format(self.C_rad.dtype.kind, sp.int_(1).dtype.kind)
        assert s_row_range.shape == (shape[0] - 1 + 2 * self.S_rad + 2 * self.C_rad + 1, ), '''s_row_range: {0}
rl: {1}; rh: {2}
wanted: {3}'''.format(s_row_range, rl, rh, (shape[0] - 1 + 2 * self.S_rad + 2 * self.C_rad + 1, ))
        assert s_col_range.shape == (shape[1] - 1 + 2 * self.S_rad + 2 * self.C_rad + 1, ), '''s_col_range: {0}
cl: {1}; ch: {2}
wanted: {3}'''.format(s_col_range, cl, ch, (shape[1] - 1 + 2 * self.S_rad + 2 * self.C_rad + 1, ))
        assert rc_row_range.shape == (shape[0] - 1 + 2 * self.C_rad + 1, ), '''rc_row_range: {0}
rl: {1}; rh: {2}
wanted: {3}'''.format(rc_row_range, rl, rh, (shape[0] - 1 + 2 * self.C_rad + 1, ))
        assert rc_col_range.shape == (shape[1] - 1 + 2 * self.C_rad + 1, ), '''rc_col_range: {0}
cl: {1}; ch: {2}
wanted: {3}'''.format(rc_col_range, cl, ch, (shape[1] - 1 + 2 * self.C_rad + 1, ))

        R_sub = self.board[rc_row_range, :, 0][:, rc_col_range]
        S_sub = self.board[s_row_range, :, 1][:, s_col_range]
        C_sub = self.board[rc_row_range, :, 2][:, rc_col_range]

        assert R_sub.shape == tuple(sp.array(shape)+2), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+2)
        assert S_sub.shape == tuple(sp.array(shape)+4), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+4)
        assert C_sub.shape == tuple(sp.array(shape)+2), 'R_sub.shape: {0}\nshape: {1} and shape+2: {2}'.format(R_sub.shape, shape, sp.array(shape)+2)

        # we count how many signallers are within each cell's neighbourhood
        self.count_neighbors(s_count,
                             p1[0],
                             s_rows,
                             s_cols,
                             gene=0,
                             allele=1,
                             radius=self.S_rad)

        # a cell will produce common goods if it's receptive and cooperative and signal in its neighborhood is above threshold
        # or when it's unreceptive and cooperative, with no regard to signal in its neighbourhood.
        cooping_cells = ((C_sub & R_sub) & (s_count >= self.S_th)) | (C_sub & (R_sub ^ True))

        assert_ndim(cooping_cells, 2)
        assert (cooping_cells.shape == (3, 4) and shape == (1, 2)) or \
               (cooping_cells.shape == (4, 3) and shape == (2, 1)), \
            '''cooping_cells.shape: {0} shape: {1}'''.format(cooping_cells.shape, shape)

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
        twocellpos_r = sp.arange(rl, rh + 1) % self.board_size
        twocellpos_c = sp.arange(cl, ch + 1) % self.board_size

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
        M += (G ^ True) * Total_cost_board
        assert_ndim(M, 2)
        M = self.metabolic_baseline / M
        assert_ndim(M, 2)
        score1 = p1 * M.item(0)  # score1 is the first position's score
        score2 = p2 * M.item(1)  # score2 is the second position's score
        if shape == (2, 1):
            assert M.shape == (2, 1), 'M.shape == {0}\nWanted == {1}'.format(M.shape, (2, 1))
            if c_pos_2[0] > c_pos_1[0]:
                # their position is like this:
                # 2
                # 1
                if score1 > score2:
                    #print("comp 2 wins")
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
                else:
                    #print("comp 1 wins")
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1
                # 2
                if score1 > score2:
                    #print("comp 1 wins")
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
                else:
                    #print("comp 2 wins")
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
        else:
            assert M.shape == (1, 2), 'M.shape == {0}\nWanted == {1}'.format(M.shape, (1, 2))
            if c_pos_2[1] < c_pos_1[1]:
                # their position is like this:
                # 2 1
                if score1 > score2:
                    #print("comp 2 wins")
                    # competitor 2 wins
                    return (c_pos_2t, c_pos_1)
                else:
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
            else:
                # their position is like this:
                # 1 2
                if score1 > score2:
                    #print("comp 1 wins")
                    # competitor 1 wins
                    return (c_pos_1, c_pos_2t)
                else:
                    #print("comp 2 wins")
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

    def save_h5(self):
        """
        Saves the attributes of self whose names show up as keys in self.parameters.
        """
        with h5py.File(self.data_filename) as ff:
            for key in self.parameters:
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
