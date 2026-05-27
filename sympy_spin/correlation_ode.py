from sympy_spin.neel import MzOperator, NzOperator, CpOperator, CmOperator, NpOperator, NmOperator, NeelOperator
from sympy_spin.spin import Sz
from functools import reduce
import sympy as sym
from sympy.physics.quantum import Commutator, Operator
from sympy_spin.common import hbar
from sympy_spin.brackets import SymmetricBracket


def decouple2(expr):
    if isinstance(expr, (sym.Number, sym.Symbol, Operator, sym.Function)) and not isinstance(expr, NeelOperator):
        return expr

    if isinstance(N := expr, NeelOperator):
        stem = N.name.split("}^")[0].split("hat{")[1]
        component = N.name.split("^")[1].split("_")[0]
        bond = N.name.split("_{")[1].split("}}")[0]

        if component == "+":
            component = "p"

        if component == "-":
            component = "m"

        i, j = bond.split(",")
        label = f"{stem}{component}_{i}_{j}"
        return sym.Symbol(label)

    else:
        return expr.func(*(decouple2(node) for node in expr.args))


def get_decoupling_map(lattice, decoupling_pattern=decouple2):
    dec_map = {}

    for O in [NzOperator, MzOperator]:
        for bond in lattice.nearest_bonds + lattice.pbc_bonds:
            dec_map[O(*bond)] = decoupling_pattern(O(*bond))

    for O in [NpOperator, NmOperator]:
        for bond in lattice.inter_bonds:
            dec_map[O(*bond)] = decoupling_pattern(O(*bond))

    for O in [CpOperator, CmOperator]:
        for dbond in lattice.intra_bonds:
            dec_map[O(*dbond)] = decoupling_pattern(O(*dbond))

    return dec_map


def get_Sz_substitution_map_uniform(lattice):
    """Return a map to replace all Sz uniformly by the Neel operators on all connected bonds."""
    bonds = lattice.nearest_bonds + lattice.pbc_bonds

    sz_map = {}

    for j in lattice.indices_a:
        bonds_containing_n = [(k, l) for (k, l) in bonds if j in [k, l]]

        sum_over_bonds = reduce(sym.Add, [MzOperator(*bond) + NzOperator(*bond) for bond in bonds_containing_n])
        Sz_shared = (sum_over_bonds / len(bonds_containing_n)).simplify()
        sz_map[Sz(j)] = Sz_shared

    for j in lattice.indices_b:
        bonds_containing_n = [(k, l) for (k, l) in bonds if j in [k, l]]

        sum_over_bonds = reduce(sym.Add, [MzOperator(*bond) - NzOperator(*bond) for bond in bonds_containing_n])
        Sz_shared = (sum_over_bonds / len(bonds_containing_n)).simplify()
        sz_map[Sz(j)] = Sz_shared
    return sz_map


def get_rhs_expressions(hamiltonian, lattice, neel_container=None, verbose=False, progress=True, eta=0, check_correctness=False):
    rhs_expressions = []

    sz_replacements = get_Sz_substitution_map_uniform(lattice)
    dec_replacements = get_decoupling_map(lattice)

    if neel_container:
        dec_replacements_ode = {O: neel_container[idx] for idx, O in enumerate(dec_replacements)}
        dec_replacements = dec_replacements_ode

    cinverted_replacements = dict(
        [(CpOperator(*reversed(dbond)), CmOperator(*dbond)) for dbond in lattice.intra_bonds] + \
        [(CmOperator(*reversed(dbond)), CpOperator(*dbond)) for dbond in lattice.intra_bonds])

    num_equations = len(dec_replacements)

    bond_classes = [
        lattice.nearest_bonds + lattice.pbc_bonds,
        lattice.inter_bonds,
        lattice.intra_bonds
    ]

    operators_for_bond_classes = [
        [NzOperator, MzOperator],
        [NpOperator, NmOperator],
        [CpOperator, CmOperator],
    ]

    progress = 0

    if check_correctness:
        H_spin = hamiltonian.subs_neel()

    for bond_class, operators in zip(bond_classes, operators_for_bond_classes):
        for O in operators:
            for bond in bond_class:

                if verbose:
                    print(f"{O.__name__}{bond}")

                dNdt = (sym.I / hbar * Commutator(hamiltonian, O(*bond)).expand(commutator=True).expand(
                    commutator=True).expand(commutator=True).expand(commutator=True).expand(commutator=True).expand(
                    commutator=True).expand(commutator=True).expand(commutator=True).apply_commutators())

                if eta != 0:
                    dNdt -= eta * SymmetricBracket(hamiltonian, O(*bond)).expand(commutator=True).expand(
                        commutator=True).expand(commutator=True).expand(commutator=True).expand(commutator=True).expand(
                        commutator=True).apply_symmetric_brackets().simplify()

                dNdt = dNdt.xreplace(sz_replacements | cinverted_replacements)

                rhs = dNdt  # .subs(sz_replacements).expand()

                if check_correctness:
                    dNdt_spin = (I / hbar * Commutator(H_spin, O(*bond).subs_neel()).expand(commutator=True).expand(
                        commutator=True).expand(commutator=True).expand(commutator=True).expand(commutator=True).expand(
                        commutator=True).expand(commutator=True).expand(commutator=True).apply_commutators())

                    difference = (dNdt_spin.expand().order_spins() - rhs.subs_neel().expand().order_spins()).simplify()
                    difference = difference.subs(
                        [(Sz(i) ** 2, S * (S + 1) * hbar ** 2 - Sy(i) ** 2 - Sx(i) ** 2) for i in
                         range(lattice.nx * lattice.ny)]).simplify()
                    difference = difference.expand().order_spins().simplify()

                    if difference != 0:
                        display(difference)

                dec_replacements[hbar] = 1
                rhs = rhs.xreplace(dec_replacements)  # .expand()
                rhs_expressions.append(rhs)

                if verbose:
                    display(rhs)

                progress += 1

                output_str = f"progress: {progress}/{num_equations}, correct: {difference == 0}" if check_correctness else f"progress: {progress}/{num_equations}"
                print(output_str, end="\r")
    return rhs_expressions
