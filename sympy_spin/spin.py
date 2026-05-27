from sympy.physics.quantum import Operator,Dagger
from sympy import I, LeviCivita, KroneckerDelta, Matrix, sqrt

from sympy_spin.common import COMPONENT_NUMBERS, hbar, Component, S
from sympy_spin.boson import BosonOperator

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

    def symmetric_bracket_with(self, other_spin):
        bracket = - KroneckerDelta(self.index, other_spin.index) * self * other_spin

        if self.component == other_spin.component:
            bracket += KroneckerDelta(self.index, other_spin.index) * S*(S+1)*hbar**2

        return bracket

    def bosonize_to_HP(self, linear):
        a = BosonOperator(self.index)

        if self.component is Component.Z:
            return hbar*(S-Dagger(a)*a)

        Sp = hbar*sqrt(2*S)*a if linear else sqrt(2*S-Dagger(a)*a)*a
        Sm = hbar*sqrt(2*S)*Dagger(a) if linear else Dagger(a)*sqrt(2*S-Dagger(a)*a)

        if self.component is Component.X:
            return ((Sp + Sm)/2).simplify()

        if self.component is Component.Y:
            return ((Sp - Sm) / 2*I).simplify()

    def has_symmetric_bracket_with(self, operator):
        return isinstance(operator, SpinOperator)

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


### Matrices in Sz basis

def Sx_matrix(S):
    
    match S:
        case 0.5:
            return hbar/2 * Matrix([[0, 1], [1, 0]])
        
        case 1:
            return hbar/sqrt(2)  * Matrix([
                                           [0, 1, 0], 
                                           [1, 0, 1],
                                           [0, 1, 0]
                                           ])

        case 1.5:
            return hbar/(2)  * Matrix([
                                           [0, sqrt(3), 0, 0], 
                                           [sqrt(3), 0, 2, 0],
                                           [0, 2, 0, sqrt(3)],
                                           [0, 0, sqrt(3), 0]
                                           ])
        
        case 2:
            return hbar/(2)  * Matrix([
                                           [0, 2, 0, 0, 0], 
                                           [2, 0, sqrt(6), 0, 0],
                                           [0, sqrt(6), 0, sqrt(6), 0],
                                           [0, 0, sqrt(6), 0, 2],
                                           [0, 0, 0, 2, 0]
                                           ])
        
        

def Sy_matrix(S):
    
    match S:
        case 0.5:
            return hbar/(2*I) * Matrix([[0, 1], [-1, 0]])
        
        case 1:
            return hbar/(sqrt(2)*I)  * Matrix([
                                           [0, 1, 0], 
                                           [-1, 0, 1],
                                           [0, -1, 0]
                                           ])
        
        case 1.5:
            return hbar/(2*I)  * Matrix([
                                           [0, sqrt(3), 0, 0], 
                                           [-sqrt(3), 0, 2, 0],
                                           [0, -2, 0, sqrt(3)],
                                           [0, 0, -sqrt(3), 0]
                                           ])
        case 2:
            return hbar/(2*I)  * Matrix([
                                           [0, 2, 0, 0, 0], 
                                           [-2, 0, sqrt(6), 0, 0],
                                           [0, -sqrt(6), 0, sqrt(6), 0],
                                           [0, 0, -sqrt(6), 0, 2],
                                           [0, 0, 0, -2, 0]
                                           ])
        

def Sz_matrix(S):
    
    match S:
        case 0.5:
            return hbar/2 * Matrix([[1, 0], [0, -1]])
        
        case 1:
            return hbar  * Matrix([
                                [1, 0, 0], 
                                [0, 0, 0],
                                [0, 0, -1]
                                ])
        
        case 1.5:
            return hbar/(2)  * Matrix([
                                           [3, 0, 0, 0], 
                                           [0, 1, 0, 0],
                                           [0, 0, -1, 0],
                                           [0, 0, 0, -3]
                                           ])
        case 2:
            return hbar  * Matrix([
                                           [2, 0, 0, 0, 0], 
                                           [0, 1, 0, 0, 0],
                                           [0, 0, 0, 0, 0],
                                           [0, 0, 0, -1, 0],
                                           [0, 0, 0, 0, -2]
                                           ])

def Sp_matrix(S):
    
    match S:
        case 0.5:
            return hbar * Matrix([[0, 1], [0, 0]])
        
        case 1:
            return hbar*sqrt(2)  * Matrix([
                                           [0, 1, 0], 
                                           [0, 0, 1],
                                           [0, 0, 0]
                                           ])
        case 1.5:
            return hbar  * Matrix([
                                           [0, sqrt(3), 0, 0], 
                                           [0, 0, 2, 0],
                                           [0, 0, 0, sqrt(3)],
                                           [0, 0, 0, 0]
                                           ])
        case 2:
            return hbar  * Matrix([
                                           [0, 2, 0, 0, 0], 
                                           [0, 0, sqrt(6), 0, 0],
                                           [0, 0, 0, sqrt(6), 0],
                                           [0, 0, 0, 0, 2],
                                           [0, 0, 0, 0, 0]
                                           ])

        
def Sm_matrix(S):
    
    match S:
        case 0.5:
            return hbar * Matrix([[0, 0], [1, 0]])
        case 1:
            return hbar*sqrt(2)  * Matrix([
                                           [0, 0, 0], 
                                           [1, 0, 0],
                                           [0, 1, 0]
                                           ])
        case 1.5:
            return hbar  * Matrix([
                                           [0, 0, 0, 0], 
                                           [sqrt(3), 0, 0, 0],
                                           [0, 2, 0, 0],
                                           [0, 0, sqrt(3), 0]
                                           ])
        
        case 2:
            return hbar  * Matrix([
                                           [0, 0, 0, 0, 0], 
                                           [2, 0, 0, 0, 0],
                                           [0, sqrt(6), 0, 0, 0],
                                           [0, 0, sqrt(6), 0, 0],
                                           [0, 0, 0, 2, 0]
                                           ])