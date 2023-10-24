import numpy as np
import matplotlib.pyplot as plt


filenames = ['data_folder/fcc100a108.txt', 'data_folder/fcc100a256.txt', 'data_folder/fcc100a500.txt', 'data_folder/fcc100a864.txt', 'data_folder/fcc100a1372.txt', 'data_folder/fcc100a2048.txt']
positions = [np.genfromtxt(filenames[i]) for i in range(len(filenames))]

def relative_pos(pos):
    n = int(np.size(pos, 0))
    rij = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            rij[i][j] = np.sqrt(np.sum((pos[i] - pos[j])**2))
    return rij

def lenny_j(r, sigma, epsilon): #sigma = 2.644 Å, epsilon = 0.345 eV
    LJ = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return LJ


def potential_energy(pos, sigma, epsilon):
    r_ij = relative_pos(pos)
    rows, cols = r_ij.shape
    s = sigma
    e = epsilon
    V = 0
    for i in range(rows):
        for j in range(i+1, cols):
            V += lenny_j(r_ij[i, j], s, e)
    return V

V_values = np.zeros((len(filenames), 1))
for k in range(len(filenames)):
    rows, cols = positions[k].shape
    sigma = 2.644
    epsilon = 0.345 
    V_f = potential_energy(positions[k], sigma, epsilon)/rows
    V_values[k] = V_f


print(V_values)
x_Vals = [108, 256, 500, 864, 1372, 2048]

fig, ax = plt.subplots()
ax.plot(x_Vals, V_values)
plt.show()
    
    
    



    
