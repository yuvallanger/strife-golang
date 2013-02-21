#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
    plt.clf()
    plt.hold(True)
    for strain in range(8):
        print 'plotting strain:', legend[strain]
        plt.plot(frequencies_samples['generation'], frequencies_samples['data'][:,strain])
    plt.ylabel('# cells of strain')
    plt.xlabel('generation #')
    plt.legend(legend)
    output_filename = output_filename_prefix+'-frequencies.png'
    print 'Saving figure:', output_filename
    plt.savefig(output_filename)
    plt.hold(False)

def get_allele_frequencies(frequencies_samples):
    r_freqs = frequencies_samples['data'][:,1] + \
              frequencies_samples['data'][:,3] + \
              frequencies_samples['data'][:,5] + \
              frequencies_samples['data'][:,7]
               
    s_freqs = frequencies_samples['data'][:,2] + \
              frequencies_samples['data'][:,3] + \
              frequencies_samples['data'][:,6] + \
              frequencies_samples['data'][:,7]

    g_freqs = frequencies_samples['data'][:,4] + \
              frequencies_samples['data'][:,5] + \
              frequencies_samples['data'][:,6] + \
              frequencies_samples['data'][:,7]

    allele_freqs = np.empty((frequencies_samples['data'].shape[0], 3), dtype='int64')

    allele_freqs[:,0] = r_freqs
    allele_freqs[:,1] = s_freqs
    allele_freqs[:,2] = g_freqs

    return allele_freqs

def output_allele_frequencies(output_filename_prefix, frequencies_samples):
    print "Starting working on allele frequencies"
    allele_freqs = get_allele_frequencies(frequencies_samples)
    plt.hold(False)
    plt.clf()
    plt.hold(True)
    
    for gene in range(3):
        print "Plotting", gene_legend[gene]
        plt.plot(frequencies_samples['generation'], allele_freqs[:,gene])
    plt.ylabel('# cells with functioning')
    plt.xlabel('generation #')
    plt.legend(gene_legend)
    output_filename = output_filename_prefix+'-allele-frequencies.png'
    print 'Saving image', output_filename
    plt.savefig(output_filename)
    plt.hold(False)
    plt.clf()

gene_legend = ['R', 'S', 'C']

def output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples):
    print 'Starting working on neighbors frequencies samples'
    plt.clf()
    for center_strain in range(8):
        print 'Plotting neighbors of strain', legend[center_strain]
        plt.hold(True)
        for neighbor_strain in range(8):
            print 'Plotting neighbor strain', legend[neighbor_strain]
            if center_strain != neighbor_strain:
                plt.plot(neighbors_frequencies_samples['generation'],
                            neighbors_frequencies_samples['data'][:,center_strain,neighbor_strain])
        plt.ylabel('# cells neighboring strain ' + legend[center_strain])
        plt.xlabel('generation #')
        plt.legend(legend)
        output_filename = output_filename_prefix+'-'+str(center_strain)+'-neighbors-frequencies.png'
        print 'Saving image', output_filename
        plt.savefig(output_filename)
        plt.hold(False)
        plt.clf()

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
    if frequencies_samples is not None:
        output_frequencies_samples(output_filename_prefix, frequencies_samples)
        output_allele_frequencies(output_filename_prefix, frequencies_samples)
    if neighbors_frequencies_samples is not None: output_neighbors_frequencies_samples(output_filename_prefix, neighbors_frequencies_samples)
    for n in sys.argv:
        if n == '--images':
            board_snapshot_samples = numpyify_sample_sequence(data['DataSamples']['Snapshots'])
            if board_snapshot_samples is not None: output_board_samples(output_filename_prefix, board_snapshot_samples)

if __name__ == '__main__':
    main()
