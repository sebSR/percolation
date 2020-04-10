try:
    import numpy as np                  # adding support for large, multi-dimensional arrays and matrices
except:
    lib_2 = 'Numpy'
    print('Program requires ' + lib_2 + ' python module to be installed')
    print('Recommended' + lib_2 + '1.17.4')
    exit(1)
try:
    import matplotlib.pyplot as plt     # Python 2D plotting library
    import matplotlib
except:
    lib_3 = 'Matplotlib'
    print('Program requires ' + lib_3 + ' python module to be installed')
    exit(1)
import random
import os
import time


"""class definition for perlocation problem"""
class Forest():                                       # class for perlocation model = square L by L
    def __init__(self,size=10, p=60, steps=100):
        self.N = steps                                # number of monte carlo steps to simulation
        self.n = 0                                    # counter, current number of monte carlo steps
        self.l = 0                                    # counter, l+=1, if there is a connection between two sides
        self.d_min = 0                                # minimal path of connection from burning method
        self.d_min_average = 0                        # minimal average path
        self.perlocation_p = 0                        # is equal l/n
        self.density = p                              # denisty of trees
        self.size = size                              # size of square
        self.M_max_average = 0
        self.M_average = 0
        self.grid = np.zeros((self.size, self.size))  # grid to plot of results
        self.M = np.zeros(self.size*self.size)        # the mass of cluster k, counts the number of sites belonging to a given cluster k
        self.M_M = []                                 # M_M = copy of M without zeros for calculs
        self.distribution = np.zeros(self.size*self.size)
        self.monte_carlo()


    def next(self,i):
        if i < self.size-1:
            return i+1
        else:
            return i


    def previous(self,i):
        if i > 0:
            return i-1
        else:
            return i


    def top_neighbour(self,i,j):
        if i > 0:
            return int(self.grid[i-1][j])
        else:
            return 0


    def left_neighbour(self,i,j):
        if j > 0:
            return int(self.grid[i][j-1])
        else:
            return 0


    def monte_carlo(self):                          # the most natural simulation method
        if self.density == 0.592746:                # every simulation has step 0.01
            for i in range(self.size):              # so in my impemenation 1% only for this case we make more digits in randoms
                for j in range(self.size):
                    r = random.uniform(0,100)      # the new random configuration
                    if r < self.density:
                        self.grid[i][j] = 1
                    else:
                        self.grid[i][j] = 0
        else:
            for i in range(self.size):
                for j in range(self.size):
                    r = random.randint(0,100)      # the new random configuration
                    if r < self.density:
                        self.grid[i][j] = 1
                    else:
                        self.grid[i][j] = 0


    def burning_method(self):
            t = 2                                   # initial paremeter for the burning method
            for i in range(self.size):              # set the fire
                if self.grid[i][0] == 1:
                    self.grid[i][0] = t             # fire all trees from first row
            state = True
            while state:
                burned = False
                for i in range(self.size):          # check all cells
                    for j in range(self.size):
                        if self.grid[i][j] == t:
                            if j == self.size-1:    # return True if the second side is reached
                                self.d_min = t-1    # set a minimal path between two sides
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
                if burned == False:		            # if there is no tree to burn then exit
                    return False
                t+=1


    def simulation_burning(self):
        self.l = 0
        while self.n < self.N:
            self.d_min = 0
            self.monte_carlo()                      	# new configuration for system
            if self.burning_method() == True:
                self.l+=1
            self.n+=1
            self.d_min_average+=self.d_min         	    # 0 or min path from burning method
        self.d_min_average = int(self.d_min_average/self.N)
        self.perlocation_p = self.l/self.N


    def plot(self, tekst):                              # make plot of the current state
        colors = ['#FFFFFF','#000000']            # black, white
        tmp = matplotlib.colors.ListedColormap(colors)
        plt.imshow(self.grid, cmap=tmp);
        plt.colorbar()
        plt.savefig(tekst+'.png')
        plt.close()


    def hoshen_kopelman(self):                          # Hoshen–Kopelman algorithm
        k = 2                                           # initial parameter
        for i in range(self.size):                      # loop after all sites
            for j in range(self.size):
                if self.grid[i][j] == 1:                # if site is occupied
                    top = self.top_neighbour(i,j)
                    left = self.left_neighbour(i,j)
                    if top == 0 and left == 0:          # if both neighbours are empty
                        self.grid[i][j] = k             # then create new cluster
                        self.M[k] = 1                   # and mass of new cluster is equal 1
                        k+=1
                    elif left != 0 and top == 0:        # if one of them is occupied
                        self.grid[i][j] = left
                        self.M[left]+=1
                    elif left == 0 and top != 0:
                        self.grid[i][j] = top
                        self.M[top]+=1
                    elif left != 0 and top != 0 and left == top:
                        self.grid[i][j] = left
                        self.M[left]+=1
                    else:                               # if both of them are occupied but with different k
                        self.grid[i][j] = left
                        self.M[left] = self.M[left] + self.M[top] + 1
                        self.M[top] = 0


    def simulation_hoshen_kopelman(self):               # Hoshen–Kopelman algorithm for Monte Carlo
        while self.n < self.N:
            self.monte_carlo()
            for i in range(len(self.M)):
                self.M[i] = 0
            self.hoshen_kopelman()
            self.M_M = self.M[self.M!=0]                # delete all zeros
            if len(self.M_M) != 0:
                self.M_average+=sum(self.M_M)/len(self.M_M)
                self.M_max_average+=max(self.M_M)
            self.n+=1
        self.M_max_average = self.M_max_average/self.N
        self.M_average = self.M_average/self.N


    def simulation_distribution(self):
        while self.n < self.N:
            self.monte_carlo()
            for i in range(len(self.M)):
                self.M[i] = 0                           # make zeros table of clusters
            self.hoshen_kopelman()                      # -> new configuration of clusters given by M
            for i in range(len(self.M)):
                if self.M[i] != 0:
                    self.M_M.append(int(self.M[i]))     # M_M = copy of M without zeros
            if len(self.M_M) != 0:
                for i in range(int(max(self.M_M))):
                    self.distribution[i] += self.M_M.count(i)
            self.n+=1
        self.distribution/=self.N

"""class finish"""


"""function for plot"""

def plot_perlocation_threshold():
    data1 = np.loadtxt('Ave_L_10_T_10000.dat')
    data2 = np.loadtxt('Ave_L_50_T_1000.dat')
    data3 = np.loadtxt('Ave_L_100_T_1000.dat')
    axes = plt.gca()
    axes.set_xlim([0,1])
    axes.set_ylim([-0.05,1.05])
    plt.scatter(data1[:,0],data1[:,1], s=20, marker='s',facecolors='none', edgecolors='r',label="L=10,T=10^4")
    plt.scatter(data2[:,0],data2[:,1], s=20, marker='o',facecolors='none', edgecolors='k',label="L=50,T=10^3")
    plt.scatter(data3[:,0],data3[:,1], s=20, marker='^',facecolors='none', edgecolors='b',label="L=100,T=10^3")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('Perlocation probability')
    plt.grid(linewidth=0.1)
    plt.savefig('plot_perlocation_threshold.eps')
    plt.close()


def plot_min_path():
    data1 = np.loadtxt('Ave_L_10_T_10000.dat')
    data2 = np.loadtxt('Ave_L_50_T_1000.dat')
    data3 = np.loadtxt('Ave_L_100_T_1000.dat')
    axes = plt.gca()
    axes.set_xlim([0,1])
    #axes.set_ylim([-0.05,1.05])
    plt.scatter(data1[:,0],data1[:,2], s=20, marker='s',facecolors='none', edgecolors='r',label="L=10,T=10^4")
    plt.scatter(data2[:,0],data2[:,2], s=20, marker='o',facecolors='none', edgecolors='k',label="L=50,T=10^3")
    plt.scatter(data3[:,0],data3[:,2], s=20, marker='^',facecolors='none', edgecolors='b',label="L=100,T=10^3")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average shortest path')
    plt.grid(linewidth=0.1)
    plt.savefig('plot_min_path.eps')
    plt.close()


def plot_max_cluster():
    data1 = np.loadtxt('Ave_2_L_10_T_10000.dat')
    data2 = np.loadtxt('Ave_2_L_50_T_1000.dat')
    data3 = np.loadtxt('Ave_2_L_100_T_1000.dat')
    axes = plt.gca()
    axes.set_xlim([0,1])
    plt.scatter(data1[:,0],data1[:,1], s=20, marker='.',facecolors='none', edgecolors='r',label="L=10,T=10^4")
    plt.scatter(data2[:,0],data2[:,1], s=20, marker='s',facecolors='none', edgecolors='k',label="L=50,T=10^3")
    plt.scatter(data3[:,0],data3[:,1], s=20, marker='^',facecolors='none', edgecolors='b',label="L=1000,T=10^3")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average size of maximum cluster')
    plt.grid(linewidth=0.1)
    plt.savefig('plot_maximum_cluster.eps')
    plt.close()


def plot_average_cluster():
    data1 = np.loadtxt('Ave_2_L_10_T_10000.dat')
    data2 = np.loadtxt('Ave_2_L_50_T_1000.dat')
    data3 = np.loadtxt('Ave_2_L_100_T_1000.dat')
    axes = plt.gca()
    axes.set_xlim([0,1])
    plt.scatter(data1[:,0],data1[:,2], s=20, marker='.',facecolors='none', edgecolors='r',label="L=10,T=10^4")
    plt.scatter(data2[:,0],data2[:,2], s=20, marker='s',facecolors='none', edgecolors='k',label="L=50,T=10^3")
    plt.scatter(data3[:,0],data3[:,2], s=20, marker='^',facecolors='none', edgecolors='b',label="L=100,T=10^3")
    plt.legend()
    plt.xlabel('p')
    plt.ylabel('average size of average cluster')
    plt.grid(linewidth=0.1)
    plt.savefig('plot_average_cluster.eps')
    plt.close()


def plot_distribution():
    p_1 = [0.2, 0.3, 0.4, 0.5]
    p_2 =  [0.592746]
    p_3 = [0.6, 0.7, 0.8]
    sign = ['s','o','v','p']
    L = 50
    T = 1000
    data_1 = np.loadtxt('Dist_p_' + str(0.2) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_2 = np.loadtxt('Dist_p_' + str(0.3) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_3 = np.loadtxt('Dist_p_' + str(0.4) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_4 = np.loadtxt('Dist_p_' + str(0.5) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_5 = np.loadtxt('Dist_p_' + str(0.592746) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_6 = np.loadtxt('Dist_p_' + str(0.6) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_7 = np.loadtxt('Dist_p_' + str(0.7) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    data_8 = np.loadtxt('Dist_p_' + str(0.8) + '_L_' + str(L) + '_T_' + str(T) + '.dat')
    plt.scatter(np.log(data_1[:,0]),np.log(data_1[:,1]), s=20, marker='s',facecolors='none', edgecolors='black',label="p=0.2")
    plt.scatter(np.log(data_2[:,0]),np.log(data_2[:,1]), s=20, marker='o',facecolors='none', edgecolors='blue',label="p=0.3")
    plt.scatter(np.log(data_3[:,0]),np.log(data_3[:,1]), s=20, marker='v',facecolors='none', edgecolors='red',label="p=0.4")
    plt.scatter(np.log(data_4[:,0]),np.log(data_4[:,1]), s=20, marker='p',facecolors='none', edgecolors='green',label="p=0.5")
    plt.scatter(np.log(data_5[:,0]),np.log(data_5[:,1]), s=20, marker='^',facecolors='none', edgecolors='yellow',label="p=0.592746")
    plt.scatter(np.log(data_6[:,0]),np.log(data_6[:,1]), s=20, marker='8',facecolors='none', edgecolors='cyan',label="p=0.6")
    plt.scatter(np.log(data_7[:,0]),np.log(data_7[:,1]), s=20, marker='*',facecolors='none', edgecolors='#4D0000',label="p=0.7")
    plt.scatter(np.log(data_8[:,0]),np.log(data_8[:,1]), s=20, marker='h',facecolors='none', edgecolors='#698B69',label="p=0.8")
    plt.xlabel('log(s)')
    plt.ylabel('log(n)')
    plt.legend()
    plt.grid(linewidth=0.1)
    plt.savefig('plot_distribution_L_' + str(L) + '.eps')
    plt.close()


    plt.scatter(np.log(data_1[:,0]),np.log(data_1[:,1]), s=20, marker='s',facecolors='none', edgecolors='r',label="p=0.2")
    plt.scatter(np.log(data_2[:,0]),np.log(data_2[:,1]), s=20, marker='o',facecolors='none', edgecolors='k',label="p=0.3")
    plt.scatter(np.log(data_3[:,0]),np.log(data_3[:,1]), s=20, marker='v',facecolors='none', edgecolors='b',label="p=0.4")
    plt.scatter(np.log(data_4[:,0]),np.log(data_4[:,1]), s=20, marker='p',facecolors='none', edgecolors='g',label="p=0.5")
    plt.xlabel('log(s)')
    plt.ylabel('log(n)')
    plt.legend()
    plt.grid(linewidth=0.1)
    plt.savefig('plot_distribution_1_L_' + str(L) + '.eps')
    plt.close()


    plt.scatter(np.log(data_5[:,0]),np.log(data_5[:,1]), s=20, marker='o',facecolors='none', edgecolors='k',label="p=0.592746")
    plt.xlabel('log(s)')
    plt.ylabel('log(n)')
    plt.legend()
    plt.grid(linewidth=0.1)
    plt.savefig('plot_distribution_2_L_' + str(L) + '.eps')
    plt.close()

    plt.scatter(np.log(data_6[:,0]),np.log(data_6[:,1]), s=20, marker='s',facecolors='none', edgecolors='r',label="p=0.6")
    plt.scatter(np.log(data_7[:,0]),np.log(data_7[:,1]), s=20, marker='o',facecolors='none', edgecolors='k',label="p=0.7")
    plt.scatter(np.log(data_8[:,0]),np.log(data_8[:,1]), s=20, marker='v',facecolors='none', edgecolors='b',label="p=0.8")
    plt.xlabel('log(s)')
    plt.ylabel('log(n)')
    plt.legend()
    plt.grid(linewidth=0.1)
    plt.grid(linewidth=0.1)
    plt.savefig('plot_distribution_3_L_' + str(L) + '.eps')
    plt.close()


def plot_burning():
    plot_perlocation_threshold()
    plot_min_path()


def plot_hoshen_kopelman():
    plot_max_cluster()
    plot_average_cluster()


def main_burning():
    data = np.loadtxt('perc_ini.txt')
    L = int(data[0])
    T = int(data[1])
    p = float(data[2]*100)
    p_k = float(data[3]*100)
    dp = float(data[4]*100)
    file = open('Ave_L_' + str(L) + '_T_' + str(T) + '.dat', 'w')
    start_time = time.time()
    while p <= p_k:
        tmp = Forest(L,p,T)
        #tmp.plot(str(p))
        tmp.simulation_burning()
        file.write(str(round(p/100,3)) + '\t' + str(round(tmp.perlocation_p,2)) + '\t' + str(tmp.d_min_average) + '\n')
        print(p)
        p+=dp
    stop_time = time.time()
    print(round(stop_time - start_time,3))
    file.close()



def main_hoshen_kopelman():
    data = np.loadtxt('perc_ini.txt')
    L = int(data[0])
    T = int(data[1])
    p = float(data[2]*100)
    p_k = float(data[3]*100)
    dp = float(data[4]*100)
    file = open('Ave_2_L_' + str(L) + '_T_' + str(T) + '.dat', 'w')
    start_time = time.time()
    while p <= p_k:
        tmp = Forest(L,p,T)
        tmp.simulation_hoshen_kopelman()
        file.write(str(round(p/100,3)) + '\t' + str(int(tmp.M_max_average)) + '\t' + str(int(tmp.M_average)) + '\n')
        print(p)
        p+=dp
    stop_time = time.time()
    print(round(stop_time - start_time,3))
    file.close()


def main_distribution_of_clusters():
    data = np.loadtxt('perc_ini.txt')
    L = int(data[0])
    T = int(data[1])
    p = [0.2, 0.3, 0.4, 0.5, 0.592746, 0.6, 0.7, 0.8]
    start_time = time.time()
    for i in p:
        file = open('Dist_p_' + str(i) + '_L_' + str(L) + '_T_' + str(T) + '.dat', 'w')
        tmp = Forest(L,i*100,T)
        tmp.simulation_distribution()
        for j in range(len(tmp.distribution)):
            if tmp.distribution[j] != 0:
                file.write(str(j) + '\t' + str(round(tmp.distribution[j],3)) + '\n')
        print(i)
    stop_time = time.time()
    print(round(stop_time - start_time,3))
    file.close()



if __name__ == '__main__':
    main_burning()
    #plot_burning()
    main_hoshen_kopelman()
    #plot_hoshen_kopelman()
    main_distribution_of_clusters()
    #plot_distribution()
