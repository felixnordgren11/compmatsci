import numpy as np



class Atom:
    def __init__(self, weight, radius, position, charge=1, index=None):
        self.weight = weight
        self.index = index
        self.radius = radius
        self.position = position
        self.charge = charge
        self.num_neighbors = 0

    def set_charge(self, charge):
        self.charge = charge

    def get_position(self):
        return self.position

    def __str__(self):
        return f"{self.element} at {self.position} with charge {self.charge}"