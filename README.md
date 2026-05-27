# SympySpin

Small package built ontop of Sympy that adds additional functionality for quantum magnetism.

Features:
- Spin operators with lattice index and automatic commutation relations
- Spin-Spin correlation operators between different lattice sites, as well as code to generate their equations of motion, as introduced in [Körber et al. (2026)](https://doi.org/10.1103/b2ny-jh15)
- Boson operators with their corresponding commutation relations
- additional expression manipulation methods such as ```.apply_commutators()``` that will insert commutators if available, or ```.order_spins()``` that will order spin operators using the commutation relations


Planned features:
- Better documentation
- Expansions for spin operators into bosons (Holstein-Primakoff, Schwinger, etc.)

## Getting started

For now, simply install the package using

```bash
pip install git+https://github.com/lkoerber/sympy_spin.git
```

This package is built on-top of ``SymPy`` and includes additionial functionality in quick-and-dirty way by simply dynamically adding methods to ``SymPy`` base classes. It is intended to work in a notebook environment. Therefore, simply start each notebook with

```python
from sympy import *
from sympy_spin import *
```

and you are good to go.

## Usage

Evaluate commutation relations directly from operators

```python
i = 0
j = 1
Sx(i).commutator_with(Sy(j))  # = delta_ij * hbar* Sz(i) 
```

or equivalently with

```python
Commutator(Sx(i), Sy(j)).apply_commutators()
```

The second way is preferred as it can be used for operator products, for example

```python
Commutator(Sx(i)*Sy(i), Sz(j)).expand(commutator=True).apply_commutators()
```

This, however requires a (sometimes repeated) call of the `.expand()` method with attribute `commutator=True`. The same 
is also possible for the correlation operators. 

In order to check operator expressions for equality, you can use, for example, the `.order_spins()` method

More examples and applications are found in [Körber et al. (2026)](https://doi.org/10.1103/b2ny-jh15) and the associated research data therein.
## Cite

If you use this package in order to calculate the dynamics of spin correlations for your research, 
please consider citing [Körber et al. (2026)](https://doi.org/10.1103/b2ny-jh15)
