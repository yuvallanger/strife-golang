#!/usr/bin/env python
import numpy as np
import pylab as pl
import json
import sys
import PIL
import PIL.Image


def numpyify_sample_sequence(sequence):
    print 'Converting board data into usable numpy arrays'
    l_generation = []
    l_data = []
    if sequence is not None:
        for i in range(len(sequence)):
            l_generation.append(sequence[i]['Generation'])
            l_data.append(sequence[i]['Data'])
        return {'generation': np.array(l_generation),
                'data': np.array(l_data)}
    return None

legend = ['rsc','Rsc','rSc','RSc','rsC','RsC','rSC','RSC']

def output_frequencies_samples(output_filename_prefix, frequencies_samples):
    print 'Starting working on frequencies samples'
    pl.plt.clf()
    pl.plt.hold(True)
    for strain in range(8):
        print 'plotting strain:', legend[strain]
        pl.plt.plot(frequencies_samples['generation'], frequencies_samples['data'][:,strain])
    pl.plt.ylabel('# cells of strain')
    pl.plt.xlabel('generation #')
    pl.plt.legend(legend)
    output_filename = output_filename_prefix+'-frequencies.png'
    print 'Saving figure:', output_filename
    pl.plt.savefig(output_filename)
    pl.plt.hold(False)

def output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples):
    print 'Starting working on neighbors frequencies samples'
    pl.plt.clf()
    for center_strain in range(8):
        print 'Plotting neighbors of strain', legend[center_strain]
        pl.plt.hold(True)
        for neighbor_strain in range(8):
            print 'Plotting neighbor strain', legend[neighbor_strain]
            if center_strain != neighbor_strain:
                pl.plt.plot(neighbors_frequencies_samples['generation'],
                            neighbors_frequencies_samples['data'][:,center_strain,neighbor_strain])
        pl.plt.ylabel('# cells neighboring strain ' + legend[center_strain])
        pl.plt.xlabel('generation #')
        pl.plt.legend(legend)
        output_filename = output_filename_prefix+'-'+str(center_strain)+'-neighbors-frequencies.png'
        print 'Saving image', output_filename
        pl.plt.savefig(output_filename)
        pl.plt.hold(False)
        pl.plt.clf()

def output_board_samples(output_filename_prefix, board_snapshot_samples):
    print 'Working on board images'
    for i in range(len(board_snapshot_samples['data'])):
        image_filename = '{}-{}-{:05d}.png'.format(output_filename_prefix,'board-snapshot',board_snapshot_samples['generation'][i])
        print "Working on", image_filename
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
    print 'Using data from:', data_filename
    data = json.load(open(data_filename))

    output_filename_prefix = 'plots-' + data_filename

    frequencies_samples = numpyify_sample_sequence(data['DataSamples']['Frequencies'])
    neighbors_frequencies_samples = numpyify_sample_sequence(data['DataSamples']['NeighborsFrequencies'])
    if frequencies_samples is not None: output_frequencies_samples(output_filename_prefix, frequencies_samples)
    if neighbors_frequencies_samples is not None: output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples)
    for n in sys.argv:
        if n == '--images':
            board_snapshot_samples = numpyify_sample_sequence(data['DataSamples']['Snapshots'])
            if board_snapshot_samples is not None: output_board_samples(output_filename_prefix, board_snapshot_samples)

if __name__ == '__main__':
    main()
