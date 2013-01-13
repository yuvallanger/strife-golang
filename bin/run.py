#!/usr/bin/env python
import pyximport; pyximport.install()
import strife.model
#import scipy as sp
import sys
import os
#import time
import strife.saveload
import numpy as np
import ConfigParser
import cProfile

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
        'board_size':        10, # The length of the board. The board will be of shape (board_size, board_size).
        'generations':       10, # Each generation involves, on average, all cells in a competition, meaning
                                 #    board_size**2 competitions.
        'S_rad':              1, # Radius of the signal's effect.
        'G_rad':              1,
        'samples_per_gen':    1,
        'initial_receptives_amount': 0.5,
        'initial_signallers_amount': 0.5,
        'initial_cooperators_amount': 0.5,
        'data_filename': 'default.npz'}

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

def go(a):
#    t = time.time()
#    every = 30*60
    # TODO: Maybe add f and d somehow like in printf? {0}f {1}d
#    print "t: {0:f}, steps thus far: {1}".format(t, a.step_count)
#    steps_a = a.step_count
    while a.step_count <= 5: #a.generations:
        print a.step_count
        a.nextstep()
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

help_string = '''
{0}: A model

{0} <-h>
    <-d/--datefile> [data_filename.h5] <-c/--config> [config_filename]
'''

default_config_filename = 'default.conf'

if __name__ == '__main__':
    data_filename = default_config['data_filename']
    conf_filename = None
    if os.path.exists(default_config_filename):
        conf_filename = default_config_filename
    for i, arg in enumerate(sys.argv):
        if arg in ('--help', '-h'):
            raise help_string
        if arg in ('--data', '-d'):
            print arg, sys.argv[i+1]
            data_filename = sys.argv[i+1]
        if arg in ('--config', '-c'):
            print arg, sys.argv[i+1]
            conf_filename = sys.argv[i+1]
    config = load_config(conf_filename, default_config)
    config['data_filename'] = data_filename
    print conf_filename, data_filename
    if os.path.exists(config['data_filename']):
        print os.path.exists(config['data_filename'])
        game = strife.model.Strife(config)
        game.set_state({key: val for key, val in np.load(config['data_filename']).items()})
        #go(game)
        cProfile.run('go(game)', 'profiling_strife.pyx.prof')
    else:
        game = strife.model.Strife(config)
        game.__init__(config)
        np.savez(config['data_filename'], **game.get_state())
        print 'go!' * 10
        #go(game)
        cProfile.run('go(game)', 'profiling_strife.pyx.prof')
#    game.save_h5()
    sys.exit(0)
    
