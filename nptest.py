import numpy as np
import matplotlib.pyplot as plt
import random


filenames = ['data_folder/fcc100a108.txt', 'data_folder/fcc100a256.txt', 'data_folder/fcc100a500.txt', 'data_folder/fcc100a864.txt', 'data_folder/fcc100a1372.txt', 'data_folder/fcc100a2048.txt']
positions = [np.genfromtxt(filenames[i]) for i in range(len(filenames))]

random.seed(-325420)
filenumber = 1
n, _ = positions[filenumber].shape
k = 1/11603
rc = 4.5
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
    n1 = int(np.size(pos, 0))
    rij = np.zeros((n1,n1))
    for i in range(n1):
        for j in range(i, n1):
            if i != j:
                dx[i][j] = pos[i][0] - pos[j][0]
                dy[i][j] = pos[i][1] - pos[j][1]
                dz[i][j] = pos[i][2] - pos[j][2]
                rij[i][j] = np.sqrt((dx[i][j])**2+(dy[i][j])**2+(dz[i][j])**2)
                dx[j][i] = -dx[i][j]
                dy[j][i] = -dy[i][j]
                dz[j][i] = -dz[i][j]
                rij[j][i] = -rij[i][j]
    return [dx, dy, dz, rij]
        
    

def lenny_j(r, sigma, epsilon): #sigma = 2.644 Å, epsilon = 0.345 eV
    LJ = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return LJ


def potential_energy(pos, sigma = sigma, epsilon = epsilon):
    r_ij = relative_pos(pos)
    rows, cols = r_ij.shape
    s = sigma
    e = epsilon
    V = 0
    for i in range(rows):
        for j in range(i+1, cols):
            V += lenny_j(r_ij[i, j], s, e)
    return V



def total_potential(sigma = sigma, epsilon = epsilon, rc = rc):
    V_values = np.zeros((len(filenames), 1))
    for k in range(len(filenames)):
        rows, _ = positions[k].shape
        V_f = potential_energy(positions[k], sigma, epsilon)/rows
        V_values[k] = V_f
    return V_values



def num_neighbors(filenum = filenumber, rc = rc):
    
    # Calculate the relative positions (pairwise distances) between atoms
    relpos = relative_pos(positions[filenum])[3]
    # Initialize an array to store the number of neighbors for each atom
    numnbrs = np.zeros(n)
    
    # For each atom, count the number of neighbors within the cutoff radius
    for i in range(n):
        for j in range(i+1,n):
            if relpos[i, j] < rc:
                numnbrs[i] += 1
                numnbrs[j] += 1
    return numnbrs

def neighbors_list(filenum = filenumber, rc = rc):
    '''
    Creates a list with the neighbors for each atom in the specified cutoff radius
    '''
    # Get the number of atoms in the file
    rows, _ = positions[filenum].shape
    
    # Calculate the relative positions (pairwise distances) between atoms
    relpos = relative_pos(positions[filenum])[3]
    
    # Initialize a list to store the neighbors for each atom
    neighbors = [[] for _ in range(rows)]
    
    # For each atom, identify the neighbors within the cutoff radius
    for i in range(rows):
        for j in range(i+1, rows):
            if relpos[i][j] < rc:
                neighbors[i].append(j)
                neighbors[j].append(i)
    return neighbors

def print_atom_neighbors(atom_index, rc = rc, filenum = filenumber):
    '''
    Prints the number of neighbors for the specified atom as well as which atoms are its neighbors
    '''
    neighbors = neighbors_list(filenum, rc)
    nlist = num_neighbors(filenum, rc)
    print(f"Atom number {atom_index} has {int(nlist[atom_index])} neighbors which are: {neighbors[atom_index]} (double nmbr check: {len(neighbors[atom_index])})")

def calc_nghbrs(atom_index, rc = rc, filenum = filenumber):
    '''
    Get a tuple with number of neighbors and list of which neighbors for the specified atom with index atom_index
    '''
    neighbors = neighbors_list(filenum, rc)
    nlist = num_neighbors(filenum, rc)
    num_nbrs = int(nlist[atom_index])
    atoms_neighbors = neighbors[atom_index]
    return [num_nbrs, atoms_neighbors]

def nghb_print(num = filenumber):
    nghbrs = num_neighbors(num)
    for i in range(len(nghbrs)):
        print(f'Atom number {i} has {nghbrs[i]} neighbors')
        
#def calc_nghbrs(filenum = filenumber, rc = rc):

def ass_velocities(T_initial = 300, mass = 108*(1.66/16)*1e-27, filenum = filenumber):
    k = 1/11603 #eV/K, Boltzmann constant
    T = T_initial
    m = mass
    
    C = np.sqrt((3*k*T)/m)
    # Get the number of atoms in the file, filenum is 0-5
    n, _ = positions[filenum].shape
    velocities = np.zeros((n,3))
    for i in range(n):
        for k in range(3):
            velocities[i, k] = C * ( 2 * random.random() - 1)
    return velocities

def kinetic_energy(filenum = filenumber, mass = 108*(1.66/16)*1e-27):
    m = mass
    E_kin = 0
    velocities = ass_velocities()
    n, _ = positions[filenum].shape
    for i in range(n):
        a = 0.5 * m * ( velocities[i, 0]**2 + velocities[i, 1]**2 + velocities[i, 2]**2)
        E_kin += a
    return E_kin

def avg_temp():
    E = kinetic_energy()
    N = n
    avg_T = E*(2/3)/(N*k)
    return avg_T

def rescaled_velocities():
    _T_in = 300
    _T_act = avg_temp()
    _velocities = ass_velocities()
    _resc_velocities = np.zeros((n,3))
    
    for i in range(n):
        for j in range(3):
            _resc_velocities[i, j] = _velocities[i, j]*np.sqrt(_T_in/_T_act)
    return _resc_velocities


def calc_forces(filenum = filenumber):
    pos = positions[filenum]
    dx, dy, dz, rik = relative_pos(pos)
    f = np.zeros((n,3))
    for k in range(n):
        #print(f'index{k}')
        nghbrs_count, nghbrs_indices = calc_nghbrs(k)
        for j in range(nghbrs_count):
            i = nghbrs_indices[j]
            drik = rik[i][k]
            if rik[i][k] > 5:
                print(f'atom number {k} s, {j}th neighbor is atom number {i}')
                print(f'where they are distance {drik} from each other')
            if drik != 0:
                c1 = 24*epsilon*(sigma**6/drik**8)*((2*sigma**6)/drik**6 - 1)
                f[k,0] += c1*dx[i][k]
                f[k,1] += c1*dy[i][k]
                f[k,2] += c1*dz[i][k]
    

            

# K = kinetic_energy()
# T = avg_temp()
# a = rescaled_velocities()
# print(K, T, a)
#nghb_print()
#print_atom_neighbors(3)
a = calc_forces()
#dx, dy, dz, rik = relative_pos()
#print(rik[1])


# V_vals = total_potential()
# x_Vals = [108, 256, 500, 864, 1372, 2048]

# fig, ax = plt.subplots()
# ax.plot(x_Vals, V_vals)
# plt.show()
    
    
    



    
