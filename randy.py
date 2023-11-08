import numpy as np
import matplotlib.pyplot as plt
import random


filenames = ['data_folder/fcc100a108.txt', 'data_folder/fcc100a256.txt', 'data_folder/fcc100a500.txt', 'data_folder/fcc100a864.txt', 'data_folder/fcc100a1372.txt', 'data_folder/fcc100a2048.txt']
positions = [np.genfromtxt(filenames[i]) for i in range(len(filenames))]

random.seed(-325420)
filenumber = 1
n, _ = positions[filenumber].shape
k = 1/11603
rc = 3.0
sigma = 2.644 #Å
epsilon = 0.345 #eV

def relative_pos(pos = filenumber):
    '''
    Calculates relative positions
    '''
    pos = positions[filenumber]
    dx = np.zeros((n,n))
    dy = np.zeros((n,n))
    dz = np.zeros((n,n))
    rij = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            dx[i][j] = pos[i][0] - pos[j][0]
            dy[i][j] = pos[i][1] - pos[j][1]
            dz[i][j] = pos[i][2] - pos[j][2]
            rij[i][j] = np.sqrt((dx[i][j])**2+(dy[i][j])**2+(dz[i][j])**2)
    return rij]

a = relative_pos()
a3 = a[3]
print(a3)


