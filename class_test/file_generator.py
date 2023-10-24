import numpy as np


class SimpleCubic:
    def generate_unit_cell():
        """Generate a unit cell for a simple cubic lattice."""
        # For a simple cubic lattice, the unit cell has one atom, usually at the origin.
        return [(0, 0, 0)]

    def duplicate_unit_cell(unit_cell, size, filename):
        """Duplicate the unit cell to form the lattice and write to a file."""
        with open(filename, 'w') as f:
            for x in range(size):
                for y in range(size):
                    for z in range(size):
                        for (dx, dy, dz) in unit_cell:
                            # Calculate the position of the atom in the lattice
                            pos_x = x + dx
                            pos_y = y + dy
                            pos_z = z + dz
                            f.write(f"{pos_x} {pos_y} {pos_z}\n")

# Example usage:
sc = SimpleCubic
unit_cell = sc.generate_unit_cell()
lattice_size = 5
sc.duplicate_unit_cell(unit_cell, lattice_size, "lattice_positions.txt")