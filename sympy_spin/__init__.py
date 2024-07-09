from sympy import Basic
from sympy.physics.quantum import *

from .manipulate import order_spins, apply_commutators, order_correlators, subs_correlators
from .common import *
from .correlation_operators import *
from .spin import *


Basic.order_spins = lambda self: order_spins(self)
Basic.apply_commutators = lambda self: apply_commutators(self)
Basic.order_correlators = lambda self: order_correlators(self)
Basic.subs_correlators = lambda self: subs_correlators(self)
