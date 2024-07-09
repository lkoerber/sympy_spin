from enum import Enum

from sympy import Symbol

hbar = Symbol("hbar")
S = Symbol("S")

COMPONENT_NUMBERS = {
    "x": 1,
    "y": 2,
    "z": 3,
}


class Component(Enum):
    X = "x"
    Y = "y"
    Z = "z"

    @property
    def number(self):
        return COMPONENT_NUMBERS[self.value]
