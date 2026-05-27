from abc import abstractmethod
from sympy.physics.quantum import Operator
from sympy import I, LeviCivita, KroneckerDelta, Integer, Rational, Symbol

from sympy_spin.common import COMPONENT_NUMBERS, hbar, Component, S
from sympy_spin.spin import Sz, Sx, Sy



class NeelOperator(Operator):

    @abstractmethod
    def subs_spins(self):
        pass

    def expval(self):
        """Return the expectation value of the operator."""
        return Symbol(r"\langle" + self.name + r"\rangle")


class NzOperator(NeelOperator):

    def __init__(self, index_a, index_b):
        super().__init__()
        self.index_a = index_a
        self.index_b = index_b
        self.bond = {index_a, index_b}
        self.name = r"{\hat{N}^z_{" + str(index_a) + "," + str(index_b) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):

        shared_spins = list(self.bond & other.bond)
        if len(shared_spins) == 0:
            return Rational(0)

        if isinstance(other, (NzOperator, MzOperator)):
            return Rational(0)

        if len(shared_spins) == 1:
            if isinstance(other, NpOperator):
                return hbar / 2 * other
            if isinstance(other, NmOperator):
                return -hbar / 2 * other

            if isinstance(other, (CmOperator, CpOperator)):
                j = shared_spins[0]

                sign = +1
                if not ((j == other.index_o1 and j == self.index_a) or (j == other.index_o2 and j == self.index_b)):
                    sign *= -1

                if isinstance(other, CmOperator):
                    sign *= -1

                return sign*hbar/2*other

        if len(shared_spins) == 2:
            if isinstance(other, NpOperator):
                return hbar * other
            if isinstance(other, NmOperator):
                return -hbar * other

        raise NotImplemented("This commutator is not implemented.")
    
    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator, NmOperator))

    def symmetric_bracket_with(self, other):

        shared_spins = list(self.bond & other.bond)
        if len(shared_spins) == 0:
            return Rational(0)

        if len(shared_spins) == 1:

            if isinstance(other, NzOperator):
                return S*(S+1)*hbar**2/4 - Sz(shared_spins[0])**2/4

            if isinstance(other, MzOperator):
                sign = 1 if shared_spins[0] == self.index_a else -1
                return sign*(S*(S+1)*hbar**2/4 - Sz(shared_spins[0])**2/4)

            if isinstance(other, (NpOperator, NmOperator)):
                sign = -1 if shared_spins[0] == self.index_a else 1
                return sign*other*Sz(shared_spins[0])/2

            if isinstance(other, (CmOperator, CpOperator)):
                j = shared_spins[0]

                if self.index_a == other.index_o1:
                    return -other * Sz(j) / 2

                elif self.index_b == other.index_o2:
                    return other * Sz(j) / 2

                elif self.index_b == other.index_o1:
                    return other * Sz(j) / 2

                elif self.index_a == other.index_o2:
                    return -other * Sz(j) / 2

            else:
                raise NotImplementedError("Bracket not implemented for this operator combination.")

        if len(shared_spins) == 2:
            if isinstance(other, NzOperator):
                return -(NzOperator(self.index_a, self.index_b)**2 + MzOperator(self.index_a, self.index_b)**2)/2 + S*(S+1)*hbar**2/2

            if isinstance(other, MzOperator):
                return -self*other

            if isinstance(other, NpOperator):
                return -self*other

            if isinstance(other, NmOperator):
                return -self*other

    def subs_spins(self):
        return (Sz(self.index_a) - Sz(self.index_b)) / 2


class MzOperator(NeelOperator):

    def __init__(self, index_a, index_b):
        super().__init__()
        self.index_a = index_a
        self.index_b = index_b
        self.bond = {index_a, index_b}
        self.name = r"{\hat{M}^z_{" + str(index_a) + "," + str(index_b) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):
        shared_spins = list(self.bond & other.bond)
        if len(shared_spins) == 0:
            return Rational(0)

        if isinstance(other, (NzOperator, MzOperator)):
            return Rational(0)

        if len(shared_spins) == 1:
            sign = 1 if shared_spins[0] == self.index_a else -1

            if isinstance(other, NpOperator):
                return sign * hbar / 2 * other

            if isinstance(other, NmOperator):
                return -sign * hbar / 2 * other

            if isinstance(other, (CmOperator, CpOperator)):
                j = shared_spins[0]

                sign = +1
                if j == other.index_o2:
                    sign *= -1

                if isinstance(other, CmOperator):
                    sign *= -1

                return sign*hbar/2*other


        if len(shared_spins) == 2:
            return Rational(0)

        raise NotImplementedError("todo")

    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator, NmOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator, NmOperator))

    def symmetric_bracket_with(self, other):

        shared_spins = list(self.bond & other.bond)
        if len(shared_spins) == 0:
            return Rational(0)

        if len(shared_spins) == 1:

            if isinstance(other, MzOperator):
                return S*(S+1)*hbar**2/4 - Sz(shared_spins[0])**2/4

            if isinstance(other, NzOperator):
                return other.symmetric_bracket_with(self)

            if isinstance(other, (NmOperator, NpOperator)):
                return -other*Sz(shared_spins[0])/2


            if isinstance(other, (CmOperator, CpOperator)):
                j = shared_spins[0]

                return -other * Sz(j) / 2


            else:
                raise NotImplementedError("Bracket not implemented for this operator combination.")

        if len(shared_spins) == 2:
            if isinstance(other, NzOperator):
                return -self*other

            if isinstance(other, MzOperator):
                return -(NzOperator(self.index_a, self.index_b)**2 + MzOperator(self.index_a, self.index_b)**2)/2 + S*(S+1)*hbar**2/2

            if isinstance(other, NpOperator):
                return -self*other

            if isinstance(other, NmOperator):
                return -self*other

    def subs_spins(self):
        return (Sz(self.index_a) + Sz(self.index_b)) / 2


class NpOperator(NeelOperator):

    def __init__(self, index_a, index_b):
        super().__init__()
        self.index_a = index_a
        self.index_b = index_b
        self.bond = {index_a, index_b}
        self.name = r"{\hat{N}^+_{" + str(index_a) + "," + str(index_b) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):

        shared_spins = list(self.bond & other.bond)

        if isinstance(other, NpOperator):
            return Rational(0)

        if isinstance(other, (NzOperator, MzOperator)):
            return -other.commutator_with(self)

        if len(shared_spins) == 0:
            return Rational(0)

        if len(shared_spins) == 2:
            if isinstance(other, (NmOperator,)):
                i, j = self.index_a, self.index_b

                return hbar * NzOperator(i, j)*(S * (S + 1) * hbar ** 2 + MzOperator(i, j) ** 2 - NzOperator(i, j) ** 2)

        if len(shared_spins) == 1:

            shared_index = shared_spins[0]

            if isinstance(other, (CpOperator, CmOperator)):
                outer = [i for i in [self.index_a, self.index_b, other.index_o1, other.index_o2] if i != shared_index]

                if isinstance(other, CpOperator) and  (self.index_a == other.index_o1 or self.index_b == other.index_o2):
                    return Rational(0)

                if isinstance(other, CmOperator) and (self.index_b == other.index_o1 or self.index_a == other.index_o2):
                    return Rational(0)

                else:
                    sign = -1
                    if outer[0] == self.index_b:
                        outer = reversed(outer)
                        sign *= -1

                    return sign*hbar*NpOperator(*outer)*Sz(shared_index)


            if isinstance(other, (NmOperator,)):
                o1, o2 = other_indices = [index for index in [self.index_a, self.index_b, other.index_a, other.index_b]
                                          if index != shared_index]
                sign = 1 if shared_index == self.index_a else -1
                C_operator = CmOperator if shared_index == self.index_a else CpOperator


                return hbar*C_operator(o1, o2)*sign*Sz(shared_spins[0])




        raise NotImplementedError("todo")

    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator, NmOperator))

    def symmetric_bracket_with(self, other):

        shared_spins = list(self.bond & other.bond)

        if len(shared_spins) == 0:
            return Rational(0)

        if len(shared_spins) == 1:

            shared_index = shared_spins[0]


            if isinstance(other, (NzOperator, MzOperator)):
                return other.symmetric_bracket_with(self)

            if isinstance(other, NpOperator):
                return -self*other

            if isinstance(other, NmOperator):
                outer = [index for index in [self.index_a, self.index_b, other.index_a, other.index_b]
                         if index != shared_index]

                COperator = CpOperator if shared_index == self.index_b else CmOperator
                return S*(S+1)*hbar**2 * COperator(*outer) - self*other

            if isinstance(other, CpOperator):
                outer = [index for index in [self.index_a, self.index_b, other.index_o1, other.index_o2]
                         if index != shared_index]

                if self.index_a == other.index_o1 or self.index_b == other.index_o2:
                    return -self*other
                else:
                    if self.index_a == other.index_o2:
                        outer = reversed(outer)

                    return S * (S + 1) * hbar ** 2 * NpOperator(*outer) - self * other

            if isinstance(other, CmOperator):
                outer = [index for index in [self.index_a, self.index_b, other.index_o1, other.index_o2]
                         if index != shared_index]

                if self.index_a == other.index_o1 or self.index_b == other.index_o2:
                    if self.index_a == other.index_o1:
                        outer = reversed(outer)

                    return S*(S+1)*hbar**2 * NpOperator(*outer) - self*other
                else:
                    return -self * other

            else:
                raise NotImplementedError("Bracket not implemented for this operator combination.")

        if len(shared_spins) == 2:
            if isinstance(other, NzOperator):
                return -self*other

            if isinstance(other, MzOperator):
                return -self*other

            if isinstance(other, NpOperator):
                return -2*self**2

            if isinstance(other, NmOperator):
                return -2*self*other + S*(S+1)*hbar**2*( S*(S+1)*hbar**2-(NzOperator(self.index_a, self.index_b)**2 + MzOperator(self.index_a, self.index_b)**2))

    def subs_spins(self):
        return (Sx(self.index_a) + I*Sy(self.index_a))*(Sx(self.index_b) - I*Sy(self.index_b))/2


class NmOperator(NeelOperator):

    def __init__(self, index_a, index_b):
        super().__init__()
        self.index_a = index_a
        self.index_b = index_b
        self.bond = {index_a, index_b}
        self.name = r"{\hat{N}^-_{" + str(index_a) + "," + str(index_b) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):


        if isinstance(other, NmOperator):
            return Rational(0)

        if isinstance(other, (NzOperator, MzOperator, NpOperator)):
            return -other.commutator_with(self)

        if isinstance(other, (CpOperator, CmOperator)):
            shared_spins = list(self.bond & other.bond)

            if len(shared_spins) != 1:
                return Rational(0)

            shared_index = shared_spins[0]

            outer = [i for i in [self.index_a, self.index_b, other.index_o1, other.index_o2] if i != shared_index]

            if isinstance(other, CmOperator) and (self.index_a == other.index_o1 or self.index_b == other.index_o2):
                return Rational(0)

            if isinstance(other, CpOperator) and (self.index_b == other.index_o1 or self.index_a == other.index_o2):
                return Rational(0)

            else:
                sign = 1
                if outer[0] == self.index_b:
                    outer = reversed(outer)
                    sign *= -1

                return sign * hbar * NmOperator(*outer) * Sz(shared_index)

        raise NotImplementedError("todo")

    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NmOperator, NpOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NpOperator, NmOperator))

    def symmetric_bracket_with(self, other):
        shared_spins = list(self.bond & other.bond)

        if len(shared_spins) == 0:
            return Rational(0)

        elif len(shared_spins) == 1:

            if isinstance(other, (NzOperator, MzOperator, NpOperator)):
                return other.symmetric_bracket_with(self)

            if isinstance(other, NmOperator):
                return -self*other

            shared_index = shared_spins[0]

            if isinstance(other, CpOperator):

                outer = [index for index in [self.index_a, self.index_b, other.index_o1, other.index_o2]
                         if index != shared_index]

                if self.index_a == other.index_o1 or self.index_b == other.index_o2:
                    if self.index_a == other.index_o1:
                        outer = reversed(outer)

                    return S * (S + 1) * hbar ** 2 * NmOperator(*outer) - self * other
                else:
                    return -self * other

            if isinstance(other, CmOperator):
                outer = [index for index in [self.index_a, self.index_b, other.index_o1, other.index_o2]
                         if index != shared_index]

                if self.index_a == other.index_o1 or self.index_b == other.index_o2:
                    return -self*other
                else:
                    if self.index_a == other.index_o2:
                        outer = reversed(outer)

                    return S * (S + 1) * hbar ** 2 * NmOperator(*outer) - self * other

            else:
                raise NotImplementedError("Bracket not implemented for this operator combination.")

        elif len(shared_spins) == 2:
            if isinstance(other, NzOperator):
                return -self*other

            if isinstance(other, MzOperator):
                return -self*other

            if isinstance(other, NmOperator):
                return -2*self**2

            if isinstance(other, NpOperator):
                return -2*self*other + S*(S+1)*hbar**2*( S*(S+1)*hbar**2-(NzOperator(self.index_a, self.index_b)**2 + MzOperator(self.index_a, self.index_b)**2))

    def subs_spins(self):
        return (Sx(self.index_a) - I*Sy(self.index_a))*(Sx(self.index_b) + I*Sy(self.index_b))/2


class CpOperator(NeelOperator):

    def __init__(self, index_o1, index_o2):
        super().__init__()
        self.index_o1 = index_o1
        self.index_o2 = index_o2
        self.bond = {index_o1, index_o2}
        self.name = r"{\hat{C}^+_{" + str(index_o1) + "," + str(index_o2) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):
        if isinstance(other, (NzOperator, MzOperator, NpOperator, NmOperator)):
            return - other.commutator_with(self)

        raise NotImplementedError("todo")

    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NmOperator, NpOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NmOperator, NpOperator))

    def symmetric_bracket_with(self, other):
        if isinstance(other, (NzOperator, MzOperator, NpOperator, NmOperator)):
            return other.symmetric_bracket_with(self)

    def subs_spins(self):
        return (Sx(self.index_o1) + I*Sy(self.index_o1))*(Sx(self.index_o2) - I*Sy(self.index_o2))/2



class CmOperator(NeelOperator):

    def __init__(self, index_o1, index_o2):
        super().__init__()
        self.index_o1 = index_o1
        self.index_o2 = index_o2
        self.bond = {index_o1, index_o2}
        self.name = r"{\hat{C}^-_{" + str(index_o1) + "," + str(index_o2) + "}}"

    def _print_contents_latex(self, printer, *args):
        return self.name

    def commutator_with(self, other):
        if isinstance(other, (NzOperator, MzOperator, NpOperator, NmOperator)):
            return - other.commutator_with(self)

        raise NotImplementedError("todo")

    def has_commutator_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NmOperator, NpOperator))

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, (NzOperator, MzOperator, NmOperator, NpOperator))

    def symmetric_bracket_with(self, other):
        if isinstance(other, (NzOperator, MzOperator, NpOperator, NmOperator)):
            return other.symmetric_bracket_with(self)

    def subs_spins(self):
        return (Sx(self.index_o1) - I*Sy(self.index_o1))*(Sx(self.index_o2) + I*Sy(self.index_o2))/2


