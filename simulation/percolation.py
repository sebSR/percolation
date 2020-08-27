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


######################################################
# PERCOLATION MODEL                                  #
######################################################

"""
class definition for percolation problem: square L by L
"""
class Forest():
    def __init__(self,size=10, probability=60, trials=100):
        # number of monte carlo trials for simulation
        self.trials = trials
        # denisty of trees
        self.probability = probability
        # size of square lattice
        self.size = size
        self.current_step = 0
        self.edge_connection = 0
        # minimal path of connection (burning method)
        self.min_path = 0
        # minimal average path
        self.min_average_path = 0
        # percolation probability
        self.percolation_p = 0
        # grid to plot of results
        self.grid = np.zeros((self.size, self.size), dtype=int)
        # the mass of cluster k, counts the number of sites belonging to a given cluster k
        self.cluster_mass = np.zeros(self.size*self.size)
        self.cluster_mass_non_zeros = []
        self.average_cluster = 0
        self.max_average_cluster = 0
        self.distribution = np.zeros(self.size*self.size)
        # initial configuration
        self.monte_carlo()


    def next(self,i):
        neighbour = i+1 if (i < self.size-1) else i
        return neighbour


    def previous(self,i):
        neighbour = i-1 if i > 0 else i
        return neighbour


    def top_neighbour(self,i,j):
        neighbour = self.grid[i-1][j] if i > 0 else 0
        return neighbour


    def left_neighbour(self,i,j):
        neighbour = self.grid[i][j-1] if j > 0 else 0
        return neighbour


    # the most natural simulation method
    def monte_carlo(self):
        for i in range(self.size):
            for j in range(self.size):
                # the new random configuration
                r = random.uniform(0,100)
                type = 1 if r < self.probability else 0
                self.grid[i][j] = type


    def burning_method(self):
            # initial paremeter for the burning method
            t = 2
            # set the fire in the first line
            for i in range(self.size):
                if self.grid[i][0] == 1: self.grid[i][0] = t
            state = True
            while state:
                burned = False
                # check all cells
                for i in range(self.size):
                    for j in range(self.size):
                        if self.grid[i][j] == t:
                            if j == self.size-1:
                                self.min_path = t-1
                                return True
                            if self.grid[self.next(i)][j] == 1:
                                self.grid[self.next(i)][j] = t+1
                                burned = True
                            if self.grid[i][self.next(j)] == 1:
                                self.grid[i][self.next(j)] = t+1
                                burned = True
                            if self.grid[self.previous(i)][j] == 1:
                                self.grid[self.previous(i)][j] = t+1
                                burned = True
                            if self.grid[i][self.previous(j)] == 1:
                                self.grid[i][self.previous(j)] = t+1
                                burned = True
                 # if there is no tree to burn then exit
                if burned == False:
                    return False
                t+=1


    def simulation_burning(self):
        self.edge_connection = 0
        while self.current_step < self.trials:
            self.min_path = 0
            self.monte_carlo()
            if self.burning_method() == True:
                self.edge_connection += 1
            self.current_step += 1
            self.min_average_path += self.min_path
        self.min_average_path = self.min_average_path/self.trials
        self.percolation_p = self.edge_connection/self.trials


    # make plot of the current state
    def plot(self, tekst):
        # black & white
        colors = ['#FFFFFF','#000000']
        tmp = matplotlib.colors.ListedColormap(colors)
        plt.imshow(self.grid, cmap=tmp);
        plt.colorbar()
        plt.savefig(tekst+'.png')
        plt.close()


    # Hoshen–Kopelman algorithm
    def hoshen_kopelman(self):
        # initial parameter
        k = 2
        # loop after all sites
        for i in range(self.size):
            for j in range(self.size):
                if self.grid[i][j] == 1:
                    top = self.top_neighbour(i,j)
                    left = self.left_neighbour(i,j)
                    if top == 0 and left == 0:
                        self.grid[i][j] = k
                        self.cluster_mass[k] = 1
                        k+=1
                    elif left != 0 and top == 0:
                        self.grid[i][j] = left
                        self.cluster_mass[left]+=1
                    elif left == 0 and top != 0:
                        self.grid[i][j] = top
                        self.cluster_mass[top]+=1
                    elif left != 0 and top != 0 and left == top:
                        self.grid[i][j] = left
                        self.cluster_mass[left]+=1
                    else:
                        self.grid[i][j] = left
                        self.cluster_mass[left] = self.cluster_mass[left] + self.cluster_mass[top] + 1
                        self.cluster_mass[top] = 0


    # Hoshen–Kopelman algorithm for Monte Carlo
    def simulation_hoshen_kopelman(self):
        while self.current_step < self.trials:
            self.monte_carlo()
            for i in range(len(self.cluster_mass)):
                self.cluster_mass[i] = 0
            self.hoshen_kopelman()
            self.cluster_mass_non_zeros = self.cluster_mass[self.cluster_mass!=0]
            if len(self.cluster_mass_non_zeros) != 0:
                self.average_cluster += sum(self.cluster_mass_non_zeros)/len(self.cluster_mass_non_zeros)
                self.max_average_cluster += max(self.cluster_mass_non_zeros)
            self.current_step += 1
        self.max_average_cluster = self.max_average_cluster/self.trials
        self.average_cluster = self.average_cluster/self.trials


"""class finish"""
