import numpy as np
from atoms import Atom


class Lattice:
    
    def __init__(self, structure):
        self.atoms = []
        self.structure = structure


    def add_atom(self, atom):
        self.atoms.append(atom)

    def get_atoms(self):
        return self.atoms

    def relative_pos(self):
        pos = np.array([atom.position for atom in self.atoms])
        n = int(np.size(pos, 0))
        rij = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                rij[i][j] = np.sqrt(np.sum((pos[i] - pos[j])**2))
        return rij

    def display(self):
        for atom in self.atoms:
            print(atom)

# Example usage:
lattice_size = 3
lattice = Lattice(lattice_size, )
relative_positions = lattice.atoms
print(relative_positions)
        