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
dt = 10**(-15)
nmd = 250
m = 108*(1.66/16)*1e-27

def relative_pos(positions_file = positions[filenumber]):
    '''
    Calculates relative positions
    '''
    pos = positions_file
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
    #print([rij])
    return [dx, dy, dz, rij]

def read_positions(positions_file):
    pos = positions_file
    read_pos = np.zeros((n,3))
    #print(pos)
    for i in range(n):
        read_pos[i, 0] = pos[i, 0]
        read_pos[i, 1] = pos[i, 1]
        read_pos[i, 2] = pos[i, 2]
    return read_pos
        
        
        
####################################################
# ENERGY CALCULATIONS
####################################################



def lenny_j(r, sigma, epsilon): #sigma = 2.644 Å, epsilon = 0.345 eV
    LJ = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return LJ


def potential_energy(positions_file = positions[filenumber], sigma = sigma, epsilon = epsilon):
    r_ij = relative_pos(positions_file)[3]
    #µprint(r_ij)
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

####################################################
# NEIGHBORS CALCULATIONS
####################################################


def neighbors_list(positions_file = relative_pos(positions[filenumber]), rc = rc):
    '''
    Creates a list with the amount of neighbors for each atom in the specified cutoff radius as well as a list of the indices for each neighbor
    '''
    # Calculate the relative positions (pairwise distances) between atoms
    pos= read_positions(positions_file)
    relpos = relative_pos(pos)[3]
    
    #print(type(relpos))
    #print(relpos)
    #print(f'Relpos ij: {relpos}')
    # Initialize a list to store the neighbors for each atom
    neighbors = [[] for _ in range(n)]
    numnbrs = np.zeros(n)
    
    # For each atom, identify the neighbors within the cutoff radius
    for i in range(n):
        for j in range(i+1, n):
            
            if relpos[i][j] < rc:
                numnbrs[i] += 1
                numnbrs[j] += 1
                neighbors[i].append(j)
                neighbors[j].append(i)
    return [numnbrs, neighbors]
    nghbrs = num_neighbors(filenum)
    for i in range(len(nghbrs)):
        print(f'Atom number {i} has {nghbrs[i]} neighbors')
        

def assgn_mom_sub_velocities(T_initial = 300, mass = 108*(1.66/16)*1e-27):
    k = 1/11603 #eV/K, Boltzmann constant
    T = T_initial
    m = mass
    C = np.sqrt((3*k*T)/m)
    velocities = np.zeros((n,3))
    #assign new randomised velocities
    for i in range(n):
        for k in range(3):
            velocities[i, k] = C * ( 2 * random.random() - 1)
    avg_velocities = np.mean(velocities, axis=0)
    # Subtract the average velocities from each velocity component
    for i in range(n):
        velocities[i] -= avg_velocities
    return velocities

def kinetic_energy(velocities, mass = 108*(1.66/16)*1e-27):
    m = mass
    E_kin = 0
    for i in range(n):
        a = 0.5 * m * ( velocities[i, 0]**2 + velocities[i, 1]**2 + velocities[i, 2]**2)
        E_kin += a
    return E_kin

def avg_temp(velocities):
    E = kinetic_energy(velocities)
    N = n
    avg_T = E*(2/3)/(N*k)
    return avg_T

def rescaled_velocities(velocities):
    _T_in = 300
    _T_act = avg_temp(velocities)
    _resc_velocities = np.zeros((n,3))
    
    for i in range(n):
        for j in range(3):
            _resc_velocities[i, j] = velocities[i, j]*np.sqrt(_T_in/_T_act)
    return _resc_velocities


def calc_forces(neighbors = neighbors_list(positions[filenumber]), pos = positions[filenumber]):
    dx, dy, dz, rik = relative_pos(pos)
    f = np.zeros((n,3))
    nghbrs_count, nghbrs_indices = neighbors
    sum_c1 = 0
    
    for k in range(n):
        for j in range(int(nghbrs_count[k])):
            i = nghbrs_indices[k][j]
            drik = rik[i][k]
            if drik != 0:
                c1 = 24*epsilon*(sigma**6/drik**8)*((2*sigma**6)/drik**6 - 1)
                sum_c1 += c1
                f[k,0] += c1*dx[i][k]
                f[k,1] += c1*dy[i][k]
                f[k,2] += c1*dz[i][k]
    print(sum_c1)
    return f

def energies(positions, velocities):
    E_kin = kinetic_energy(velocities)
    E_pot = potential_energy(positions)
    E_tot = E_pot + E_kin
    return E_pot, E_kin, E_tot
    
def print_forces():
    atom_number = [i for i in range(n)]
    forces = calc_forces()
    num_neigbrs, _ = neighbors_list()
    # Print the header
    print(f"{'Atom Number /':<12}{'# of Neighbors /':<15}{'fx /':<8}{'fy /':<8}{'fz /':<8}")

    # Print data for each atom
    for i in range(n):
        force_x, force_y, force_z = forces[i]
        print(f"{atom_number[i]:<12}{num_neigbrs[i]:<15}{force_x:<8.3f}{force_y:<8.3f}{force_z:<8.3f}")


pos = positions[filenumber]
org_velocities = assgn_mom_sub_velocities()
velocities = rescaled_velocities(org_velocities)
rel_pos = relative_pos(positions[filenumber])
nghbrs_count, nghbrs_indices = neighbors_list(pos)
nghbrs = neighbors_list(pos)
ev_Energy = np.zeros((3,nmd))



for ind in range(nmd):
    forces = calc_forces(nghbrs, pos)
    rowf, colf = forces.shape
    new_forces = np.zeros((rowf,colf))
    for i in range(n):
        for j in range(3):
            pos[i][j] = pos[i][j]+velocities[i][j]*dt+0.5*(1/m)*forces[i][j]*(dt**2)
    nghbrs = neighbors_list(pos)
    new_forces = calc_forces(nghbrs, pos)
    for i in range(n):
        for j in range(3):
            velocities[i][j] = velocities[i][j]+0.5*dt*(forces[i][j]*(1/m)+new_forces[i][j])
    print(avg_temp(velocities))
    #avg_velocities = np.mean(velocities, axis=0)
    Epot, Ekin, Etot = energies(pos, velocities)
    ev_Energy[0,ind] = Epot
    ev_Energy[1,ind] = Ekin
    ev_Energy[2,ind] = Etot
    print(f'Second: {ind+1}, Avg velocity: {Ekin}')

            
    
    
        
    
#velocities[i][j] = velocities[i][j] + 0.5*dt*((1/m)*forces[i][j])


#print_forces()


# V_vals = total_potential()
# x_Vals = [108, 256, 500, 864, 1372, 2048]

# fig, ax = plt.subplots()
# ax.plot(x_Vals, V_vals)
# plt.show()
    
    
    



    
