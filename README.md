# sympy_spin
Small package built ontop of Sympy that adds additional functionality for quantum magnetism.

Features:
- Spin operators with lattice index and automatic commutation relations
- Spin-Spin correlation operators between different lattice sites
- additional expression manipulation methods such as ```.apply_commutators()``` that will insert commutators if available, or ```.order_spins()``` that will order spin operators using the commutation relations


Planned features:
- Documentation...
- Boson operators with commutation relations
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

Examples will be added here.

