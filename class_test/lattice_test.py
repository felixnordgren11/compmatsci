import numpy as np
from atoms import Atom


class Lattice:
    def __init__(self, filename):
        self.atoms = []
        self._populate_from_file(filename)

    def _populate_from_file(self, filename):
        with open(filename, 'r') as f:
            for line in f:
                x, y, z = line.split()
                atom = Atom((int(x), int(y), int(z)))
                self.add_atom(atom)

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
            
lattice = Lattice("lattice_positions.txt")
lattice.display()
rel_pos = Lattice.relative_pos()
print(rel_pos)