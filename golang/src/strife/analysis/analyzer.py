#!/usr/bin/env python
import numpy as np
import pylab as pl
import panda as pd
import json
import h5py
import sys
import PIL
import PIL.Image


def numpyify_sample_sequence(sequence):
    print 'numpyify_sample_sequence'
    print type(sequence)
    l_generation = []
    l_data = []
    if sequence is not None:
        for i in range(len(sequence)):
            l_generation.append(sequence[i]['Generation'])
            l_data.append(sequence[i]['Data'])
        return {'generation': np.array(l_generation),
                'data': np.array(l_data)}
    return None


def output_frequencies_samples(output_filename_prefix, frequencies_samples):
    print 'output_frequencies_samples'
    print type(frequencies_samples)
    pl.plt.hold(True)
    for strain in range(8):
        print 'strain:', strain, pl.plt.plot(frequencies_samples['generation'], frequencies_samples['data'][:,strain])
    pl.plt.legend(['rsc','Rsc','rSc','RSc','rsC','RsC','rSC','RSC'])
    print 'savefig:', pl.plt.savefig(output_filename_prefix+'-frequencies.png')
    pl.plt.hold(False)

def output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples):
    print 'output_neighbors_frequencies_samples'
    print type(neighbors_frequencies_samples)
    pl.plt.clf()
    for center_strain in range(8):
        pl.plt.hold(True)
        for neighbor_strain in range(8):
            print pl.plt.plot(neighbors_frequencies_samples['generation'], neighbors_frequencies_samples['data'][:,center_strain,neighbor_strain])
        pl.plt.legend(['rsc','Rsc','rSc','RSc','rsC','RsC','rSC','RSC'])
        output_filename_prefix_strain = output_filename_prefix+'-'+str(center_strain)+'-neighbors-frequencies.png'
        print output_filename_prefix_strain
        print 'savefig:', pl.plt.savefig(output_filename_prefix_strain)
        pl.plt.hold(False)
        pl.plt.clf()

def output_board_samples(output_filename_prefix, board_snapshot_samples):
    print 'output_board_samples'
    print type(board_snapshot_samples)
    for i in range(len(board_snapshot_samples['data'])):
        print i
        image_filename = '{}-{}-{:05d}.png'.format(output_filename_prefix,'board-snapshot',board_snapshot_samples['generation'][i])
        im = imagify_strain_board(board_snapshot_samples['data'][i])
        im = PIL.Image.fromarray(im, 'RGB')
        im.thumbnail((600, 600))
        im.save(image_filename)

def imagify_strain_board(strain_board):
    image_strain_board = np.empty((strain_board.shape[0], strain_board.shape[1], 3), dtype='uint8')
    for row in range(strain_board.shape[0]):
        for col in range(strain_board.shape[1]):
            image_strain_board[row, col, 0] = 255 * (strain_board[row,col] & 1) 
            image_strain_board[row, col, 1] = 255 * (strain_board[row,col] & 2)
            image_strain_board[row, col, 2] = 255 * (strain_board[row,col] & 4)
    return image_strain_board

def main():
    data_filename = sys.argv[1]
    print 'data_filename:', data_filename
    data = json.load(open(data_filename))

    output_filename_prefix = 'plots-' + data_filename

    board_snapshot_samples = numpyify_sample_sequence(data['DataSamples']['Snapshots'])
    frequencies_samples = numpyify_sample_sequence(data['DataSamples']['Frequencies'])
    neighbors_frequencies_samples = numpyify_sample_sequence(data['DataSamples']['NeighborsFrequencies'])
    print board_snapshot_samples['data'].shape
    print frequencies_samples['data'].shape
    print neighbors_frequencies_samples['data'].shape
    if board_snapshot_samples is not None: output_board_samples(output_filename_prefix, board_snapshot_samples)
    if frequencies_samples is not None: output_frequencies_samples(output_filename_prefix, frequencies_samples)
    if neighbors_frequencies_samples is not None: output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples)

if __name__ == '__main__':
    print '__name__ == "__main__"'
    main()
