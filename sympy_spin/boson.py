from sympy import KroneckerDelta, Integer
from sympy.physics.quantum import Operator


class BosonOperator(Operator):

    def __init__(self, index, species="a", is_creator=False):
        super().__init__()
        self.index = index
        self.name = r"\hat{" + str(species) + "}_{" + str(index) + "}"
        self.species = species

        self.is_creator = is_creator

    def commutator_with(self, other):
        if self.species != other.species:
            return Integer(0)

        if self.is_creator == other.is_creator:
            return Integer(0)

        sign = -1 if self.is_creator else 1
        return sign * KroneckerDelta(self.index, other.index)

    def has_commutator_with(self, operator):
        return isinstance(operator, BosonOperator)

    def _print_contents_latex(self, printer, *args):
        if self.is_creator:
            return r"{\hat{" + str(self.species) + "}_{" + str(self.index) + r"}^\dagger}"
        else:
            return r"{\hat{" + str(self.species) + "}_{" + str(self.index) + "}}"

    def _eval_adjoint(self):
        return BosonOperator(self.index, self.species, not self.is_creator)
