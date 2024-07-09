from sympy.core.operations import AssocOp
from sympy.core import Basic
from sympy.physics.quantum.commutator import Commutator, Operator
from sympy import Number, Function, Pow, Mul, Symbol

from sympy_spin.correlation_operators import CorrelationOperator
from sympy_spin.spin import SpinOperator


def expand_powers(args):
    args = list(args)
    if len(args) == 0:
        return []

    rest = [] if len(args) == 1 else args[1:]

    if isinstance(power := args[0], Pow) and not isinstance(power.base, Number):
        return [power.base for i in range(power.exp)] + expand_powers(rest)

    expanded = [args[0]] + expand_powers(rest)
    return expanded


def order_spins(expr):
    expr = expr.expand()

    def qsimplify_walktree(expr):

        if isinstance(expr, Number):
            return expr

        if not isinstance(expr, AssocOp) and not isinstance(expr, Function):
            return expr.copy()

        elif not isinstance(expr, Mul):
            return expr.func(*(qsimplify_walktree(node) for node in expr.args))

        args = expand_powers(list(expr.args))

        for i in range(len(args) - 1):

            A = args[i]
            B = args[i + 1]
            if has_commutator_with(A, B):
                if isinstance(A, SpinOperator):
                    if str(A.component) > str(B.component):
                        args = args[0:i] + [B * A + A.commutator_with(B)] + args[i + 2:]
                        result = Mul(*args).expand()
                        return qsimplify_walktree(result)
                    elif str(A.component) == str(B.component) and str(A.index) > str(B.index):
                        # print(A,B)
                        args = args[0:i] + [B * A + A.commutator_with(B)] + args[i + 2:]
                        result = Mul(*args).expand()
                        return qsimplify_walktree(result)

        return expr.copy()

    return qsimplify_walktree(expr)


def has_commutator_with(a, b):
    if not callable(getattr(a, "commutator_with", None)):
        return False

    if not callable(getattr(b, "commutator_with", None)):
        return False

    return a.has_commutator_with(b)


def apply_commutators(expr):
    if isinstance(expr, (Number, Symbol, Operator, Function)):
        return expr

    if not (func := expr.func) is Commutator:
        return expr.func(*(apply_commutators(node) for node in expr.args))

    else:
        A, B = expr.args
        if not has_commutator_with(A, B):
            return expr
        return A.commutator_with(B)


def subs_correlators(expr):
    if isinstance(expr, (Number, Symbol, Operator, Function)) and not isinstance(expr, CorrelationOperator):
        return expr

    if isinstance(C := expr, CorrelationOperator):
        return C.subs_spins()

    else:
        return expr.func(*(subs_correlators(node) for node in expr.args))


def order_correlators(expr):
    if isinstance(expr, (Number, Symbol, Operator, Function)) and not isinstance(expr, CorrelationOperator):
        return expr

    if isinstance(C := expr, CorrelationOperator):
        i = C.i
        j = C.j
        c = C.component
        if i > j:
            return -CorrelationOperator(j, i, c) if c == 2 else CorrelationOperator(j, i, c)
        return C

    else:
        return expr.func(*(order_correlators(node) for node in expr.args))