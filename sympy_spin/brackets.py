from sympy.core.add import Add
from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator

class SymmetricBracket(Expr):
    is_commutative = True

    def __new__(cls, A, B):
        r = cls.eval(A, B)
        if r is not None:
            return r
        obj = Expr.__new__(cls, A, B)
        return obj

    @classmethod
    def eval(cls, a, b):
        # if not (a and b):
        #    return S.Zero
        # if a == b:
        #    return S.Zero
        # if a.is_commutative or b.is_commutative:
        #    return S.Zero

        # [xA,yB]  ->  xy*[A,B]
        ca, nca = a.args_cnc()
        cb, ncb = b.args_cnc()
        c_part = ca + cb
        if c_part:
            return Mul(Mul(*c_part), cls(Mul._from_args(nca), Mul._from_args(ncb)))

        # Canonical ordering of arguments
        # The SymmetricBracket [A, B] is in canonical form if A < B.
        if a.compare(b) == 1:
            return cls(b, a)

    def _expand_pow(self, A, B, sign):
        exp = A.exp
        if not exp.is_integer or not exp.is_constant() or abs(exp) <= 1:
            # nothing to do
            return self
        base = A.base
        if exp.is_negative:
            base = A.base ** -1
            exp = -exp
        comm = SymmetricBracket(base, B).expand(commutator=True)

        result = base ** (exp - 1) * comm
        for i in range(1, exp):
            result += base ** (exp - 1 - i) * comm * base ** i
        return sign * result.expand()

    def _eval_expand_commutator(self, **hints):
        A = self.args[0]
        B = self.args[1]

        if isinstance(A, Add):
            # [A + B, C]  ->  [A, C] + [B, C]
            sargs = []
            for term in A.args:
                comm = SymmetricBracket(term, B)
                if isinstance(comm, SymmetricBracket):
                    comm = comm._eval_expand_commutator()
                sargs.append(comm)
            return Add(*sargs)
        elif isinstance(B, Add):
            # [A, B + C]  ->  [A, B] + [A, C]
            sargs = []
            for term in B.args:
                comm = SymmetricBracket(A, term)
                if isinstance(comm, SymmetricBracket):
                    comm = comm._eval_expand_commutator()
                sargs.append(comm)
            return Add(*sargs)
        elif isinstance(A, Mul):
            # [A*B, C] -> A*[B, C] + [A, C]*B
            a = A.args[0]
            b = Mul(*A.args[1:])
            c = B
            comm1 = SymmetricBracket(b, c)
            comm2 = SymmetricBracket(a, c)
            if isinstance(comm1, SymmetricBracket):
                comm1 = comm1._eval_expand_commutator()
            if isinstance(comm2, SymmetricBracket):
                comm2 = comm2._eval_expand_commutator()
            first = Mul(a, comm1)
            second = Mul(comm2, b)
            return Add(first, second)
        elif isinstance(B, Mul):
            # [A, B*C] -> [A, B]*C + B*[A, C]
            a = A
            b = B.args[0]
            c = Mul(*B.args[1:])
            comm1 = SymmetricBracket(a, b)
            comm2 = SymmetricBracket(a, c)
            if isinstance(comm1, SymmetricBracket):
                comm1 = comm1._eval_expand_commutator()
            if isinstance(comm2, SymmetricBracket):
                comm2 = comm2._eval_expand_commutator()
            first = Mul(comm1, c)
            second = Mul(b, comm2)
            return Add(first, second)
        elif isinstance(A, Pow):
            # [A**n, C] -> A**(n - 1)*[A, C] + A**(n - 2)*[A, C]*A + ... + [A, C]*A**(n-1)
            return self._expand_pow(A, B, 1)
        elif isinstance(B, Pow):
            # [A, C**n] -> C**(n - 1)*[C, A] + C**(n - 2)*[C, A]*C + ... + [C, A]*C**(n-1)
            return self._expand_pow(B, A, -1)

        # No changes, so return self
        return self


    def _sympyrepr(self, printer, *args):
        return "%s(%s,%s)" % (
            self.__class__.__name__, printer._print(
                self.args[0]), printer._print(self.args[1])
        )

    def _sympystr(self, printer, *args):
        return "[%s,%s]" % (
            printer._print(self.args[0]), printer._print(self.args[1]))

    def _pretty(self, printer, *args):
        pform = printer._print(self.args[0], *args)
        pform = prettyForm(*pform.right(prettyForm(',')))
        pform = prettyForm(*pform.right(printer._print(self.args[1], *args)))
        pform = prettyForm(*pform.parens(left='[', right=']'))
        return pform

    def _latex(self, printer, *args):
        return r"\left\{%s,%s\right\}" % tuple([
            printer._print(arg, *args) for arg in self.args])