try:
    # adding support for large, multi-dimensional arrays and matrices
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import random
    import os
    import time
except:
    print('Program requires python modules to be installed')
    exit(1)


"""
FUNCTION FOR PLOT
"""

def plot_percolation_threshold(trials, lattice_size):
    color_list = ['blue','red','green','yellow','cyan']
    markers_list = ['.','s','^','o','p']
    data = [np.loadtxt(f"burning_L_{i}_T_{trials}.dat") for i in lattice_size]
    axes = plt.gca()
    axes.set_xlim([0,1])
    axes.set_ylim([-0.05,1.05])
    for i in range(len(lattice_size)):
        sub_data = data[i]
        plt.scatter(sub_data[:,0],sub_data[:,1],s=20,marker=markers_list[i],facecolors='none',
                        edgecolors=color_list[i],label=f"L={lattice_size[i]},T={trials}")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('Percolation probability')
    plt.grid(linewidth=0.5)
    plt.savefig('plot_percolation_threshold.png')
    plt.close()


def plot_min_path(trials, lattice_size):
    color_list = ['blue','red','green','yellow','cyan']
    markers_list = ['.','s','^','o','p']
    data = [np.loadtxt(f'burning_L_{i}_T_{trials}.dat') for i in lattice_size]
    axes = plt.gca()
    axes.set_xlim([0,1])
    for i in range(len(lattice_size)):
        sub_data = data[i]
        plt.scatter(sub_data[:,0],sub_data[:,2], s=20, marker=markers_list[i],facecolors='none',
                    edgecolors=color_list[i],label=f"L={lattice_size[i]},T={trials}")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average shortest path')
    plt.grid(linewidth=0.5)
    plt.savefig('plot_min_path.png')
    plt.close()


def plot_max_cluster(trials, lattice_size):
    color_list = ['blue','red','green','yellow','cyan']
    markers_list = ['.','s','^','o','p']
    data = [np.loadtxt(f'clustering_L_{i}_T_{trials}.dat') for i in lattice_size]
    axes = plt.gca()
    axes.set_xlim([0,1])

    for i in range(len(lattice_size)):
        sub_data = data[i]
        plt.scatter(sub_data[:,0],sub_data[:,1], s=20, marker=markers_list[i],facecolors='none',
                    edgecolors=color_list[i],label=f"L={lattice_size[i]},T={trials}")

    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average size of maximum cluster')
    plt.grid(linewidth=0.5)
    plt.savefig('plot_maximum_cluster.png')
    plt.close()


def plot_average_cluster(trials, lattice_size):
    color_list = ['blue','red','green','yellow','cyan']
    markers_list = ['.','s','^','o','p']
    data = [np.loadtxt(f"clustering_L_{i}_T_{trials}.dat") for i in lattice_size]
    axes = plt.gca()
    axes.set_xlim([0,1])

    for i in range(len(lattice_size)):
        sub_data = data[i]
        plt.scatter(sub_data[:,0],sub_data[:,2], s=20, marker=markers_list[i],facecolors='none',
                    edgecolors=color_list[i],label=f"L={lattice_size[i]},T={trials}")

    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average size of average cluster')
    plt.grid(linewidth=0.5)
    plt.savefig('plot_average_cluster.png')
    plt.close()



def plot_burning(trials,lattice_size):
    plot_percolation_threshold(trials,lattice_size)
    plot_min_path(trials,lattice_size)


def plot_clustering(trials,lattice_size):
    plot_max_cluster(trials,lattice_size)
    plot_average_cluster(trials,lattice_size)



if __name__ == '__main__':
    # PARAMETERS FOR PLOT
    # number of steps used in simulation
    trials = 1000
    # lattice sizes which we have to done from percolation_main.py
    lattice_size = [20]

    # function plot results from burning method
    plot_burning(trials, lattice_size)

    # function plot results from Hosehen - Kopelman Algoritm
    plot_clustering(trials,lattice_size)
