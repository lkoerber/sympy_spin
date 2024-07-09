from sympy.physics.quantum import Operator
from sympy import I, LeviCivita, KroneckerDelta

from sympy_spin.common import COMPONENT_NUMBERS, hbar, Component


class SpinOperator(Operator):

    def __init__(self, index, component):
        super().__init__()
        self.index = index
        self.component = Component(str(component))
        self.name = r"\hat{S}_{" + str(self.component.value) + "," + str(index) + "}"

    def commutator_with(self, other_spin):
        comm = 0
        for c in ["x", "y", "z"]:
            number = COMPONENT_NUMBERS[c]
            comm += hbar * I * LeviCivita(self.component.number, other_spin.component.number, number) * SpinOperator(
                self.index, c) * KroneckerDelta(self.index, other_spin.index)

        # comm = Sum(I*LeviCivita(self.component, other_spin.component, k)*SpinOperator(self.index,k), (k, 1,3)).doit()
        return comm

    def has_commutator_with(self, operator):
        return isinstance(operator, SpinOperator)

    def _print_contents_latex(self, printer, *args):
        return self.name


def Sx(index):
    return SpinOperator(index, "x")


def Sy(index):
    return SpinOperator(index, "y")


def Sz(index):
    return SpinOperator(index, "z")
