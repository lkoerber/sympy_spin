from sympy import Basic
from sympy.physics.quantum import *

from .manipulate import order_spins, apply_commutators, order_correlators, subs_correlators, subs_neel, apply_symmetric_brackets,bosonize_spins_to_HP
from .common import *
from .correlation_operators import *
from .spin import *
from .brackets import *
from .boson import *
from .neel import *
from .lattices import BipartiteSquareLattice
from .correlation_ode import get_rhs_expressions, get_decoupling_map

Basic.order_spins = lambda self: order_spins(self)
Basic.apply_commutators = lambda self: apply_commutators(self)
Basic.apply_symmetric_brackets = lambda self: apply_symmetric_brackets(self)
Basic.order_correlators = lambda self: order_correlators(self)
Basic.subs_correlators = lambda self: subs_correlators(self)
Basic.subs_neel = lambda self: subs_neel(self)
Basic.bosonize_spins_to_HP = lambda self, linear: bosonize_spins_to_HP(self, linear)