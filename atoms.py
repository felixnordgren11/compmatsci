import numpy as np



class Atom:
    def __init__(self, weight, radius, position, charge=1):
        self.weight = weight
        self.radius = radius
        self.position = position
        self.charge = charge

    def set_charge(self, charge):
        self.charge = charge

    def get_position(self):
        return self.position

    def __str__(self):
        return f"{self.element} at {self.position} with charge {self.charge}"