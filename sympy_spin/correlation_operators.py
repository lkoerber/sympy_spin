from sympy.physics.quantum import Operator
from sympy import I, LeviCivita, KroneckerDelta, Integer

from sympy_spin.common import COMPONENT_NUMBERS, hbar, Component, S
from sympy_spin.spin import Sz, Sx, Sy


class CorrelationOperator(Operator):

    def __init__(self, i, j, component):

        super().__init__()
        self.i = i
        self.j = j
        self.component = component
        self.name = r"{\hat{\mathcal{C}}_{" + str(i) + str(j) + "}^{(" + str(component) + ")}}"

    def commutator_with(self, other_correlator):

        c1 = self.component
        c2 = other_correlator.component

        i = self.i
        j = self.j
        k = other_correlator.i
        l = other_correlator.j

        if ((i, j) == (k, l)) or ((i, j) == (l, k)):
            if c1 == c2:
                return Integer(0)

            if (c1, c2) in [(1, 2), (2, 1)]:
                sign = -1 if (i, j) == (k, l) else 1
                return sign * I * hbar * (
                        LeviCivita(c1, c2, 3) * (Sz(i) * C3(i, j) - C3(i, j) * Sz(j)) / (2 * S)
                ).simplify()

        if (c1, c2) == (1, 1):
            return (I * hbar / (2 * S) * (
                    (KroneckerDelta(l, i) * C2(j, k) + KroneckerDelta(i, k) * C2(j, l)) * Sz(i) + \
                    (KroneckerDelta(l, j) * C2(i, k) + KroneckerDelta(j, k) * C2(i, l)) * Sz(j)
            ))
        if (c1, c2) == (2, 2):
            return I * hbar / (2 * S) * (
                    (KroneckerDelta(i, k) * C2(j, l) - KroneckerDelta(i, l) * C2(j, k)) * Sz(i) + \
                    (-KroneckerDelta(j, k) * C2(i, l) + KroneckerDelta(j, l) * C2(i, k)) * Sz(j)
            )

        if (c1, c2) == (3, 3):
            return Integer(0)

        if (c1, c2) in [(1, 3), (3, 1)]:
            if (c1, c2) == (1, 3):
                i, j, k, l = k, l, i, j

            return I * LeviCivita(c1, c2, 2) * hbar / (2 * S) * (
                        -C2(k, l) * (Sz(j) * KroneckerDelta(i, k) + Sz(i) * KroneckerDelta(j, k)) + (
                            Sz(j) * KroneckerDelta(l, i) + Sz(i) * KroneckerDelta(l, j)) * C2(k, l))

        if (c1, c2) in [(2, 3), (3, 2)]:
            if (c1, c2) == (3, 2):
                i, j, k, l = k, l, i, j

            return I * LeviCivita(c1, c2, 1) * hbar / (2 * S) * (
                        -C1(i, j) * (Sz(l) * KroneckerDelta(i, k) + Sz(k) * KroneckerDelta(i, l)) + \
                        +(Sz(l) * KroneckerDelta(j, k) + Sz(k) * KroneckerDelta(j, l)) * C1(i, j))

        if (c1, c2) in [(1, 2), (2, 1)]:
            if (c1, c2) == (2, 1):
                i, j, k, l = k, l, i, j

            return I * LeviCivita(c1, c2, 3) * hbar / (2 * S) * (
                    Sz(i) * KroneckerDelta(l, i) * C1(j, k) + Sz(j) * KroneckerDelta(l, j) * C1(i, k) + \
                    - KroneckerDelta(i, k) * C1(j, l) * Sz(i) - KroneckerDelta(j, k) * C1(i, l) * Sz(j)
            )

        print("no commatator defined for c1 = {c1}, c2 = {c2} and i, j, k, l = {i}, {j}, {k}, {l}")
        return self * other_correlator - other_correlator * self

    def has_commutator_with(self, operator):
        return isinstance(operator, CorrelationOperator)

    def subs_spins(self):
        if self.component == 1:
            return (Sx(self.i) * Sx(self.j) + Sy(self.i) * Sy(self.j)) / (2 * S)
        if self.component == 2:
            return (Sx(self.i) * Sy(self.j) - Sx(self.j) * Sy(self.i)) / (2 * S)
        if self.component == 3:
            return Sz(self.i) * Sz(self.j) / (2 * S) + (S + 1) / 2

    def _print_contents_latex(self, printer, *args):

        return self.name


def C1(i, j):
    return CorrelationOperator(i, j, 1)


def C2(i, j):
    return CorrelationOperator(i, j, 2)


def C3(i, j):
    return CorrelationOperator(i, j, 3)
