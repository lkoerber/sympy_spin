import numpy as np
import matplotlib.pyplot as plt
from itertools import product, groupby

from io import BytesIO
import base64


COLOR_A = "#00ACB9"
COLOR_B = "#893BFF"

class BipartiteSquareLattice:
    def __init__(self, nx, ny, PBC=False):

        self.nx = nx
        self.ny = ny
        nearest_bonds = []

        grid = np.arange(0, nx * ny, 1).reshape(ny, nx)

        for vertical_bond in zip(grid[:-1].flatten(), grid[1:].flatten()):
            nearest_bonds.append(vertical_bond)

        for horizontal_bond in zip(grid[:, :-1].flatten(), grid[:, 1:].flatten()):
            nearest_bonds.append(horizontal_bond)

        self.indices = grid.flatten()

        checkerboard = np.zeros((self.ny, self.nx), dtype=bool)
        checkerboard[1::2, ::2] = True
        checkerboard[::2, 1::2] = True

        self.indices_a = grid[~checkerboard].flatten()
        self.indices_b = grid[checkerboard].flatten()

        self.nearest_bonds = []
        for i, j in nearest_bonds:
            sorted_bond = (i, j) if i in self.indices_a else (j, i)
            self.nearest_bonds.append(sorted_bond)

        self.nearest_bonds.sort()

        pbc_bonds = []
        if PBC:
            if grid.shape[0] > 2:
                pbc_bonds += list(zip(grid[0], grid[-1]))
            if len(grid.shape) >= 2 and grid.shape[1] > 2:
                pbc_bonds += list(zip(grid[:, 0], grid[:, -1]))

        self.pbc_bonds = []
        for i, j in pbc_bonds:
            sorted_bond = (i, j) if i in self.indices_a else (j, i)
            self.pbc_bonds.append(sorted_bond)

        # intralattice
        self.intra_bonds = []
        for i in self.indices_a:
            for j in self.indices_a:
                if i != j:
                    self.intra_bonds.append(tuple(sorted([i, j])))

        for i in self.indices_b:
            for j in self.indices_b:
                if i != j:
                    self.intra_bonds.append(tuple(sorted([i, j])))

        self.intra_bonds.sort()
        self.intra_bonds = list(k for k, _ in groupby(self.intra_bonds))

        self.inter_bonds = []
        inter_bonds = []
        for i in self.indices_a:
            for j in self.indices_b:
                if i != j:
                    inter_bonds.append(sorted([i, j]))

        inter_bonds.sort()
        inter_bonds = list(k for k, _ in groupby(inter_bonds))
        for i, j in inter_bonds:
            sorted_bond = (i, j) if i in self.indices_a else (j, i)
            self.inter_bonds.append(sorted_bond)

        self.not_nearest_inter_bonds = [bond for bond in self.inter_bonds if not bond in self.nearest_bonds]

        self.all_bonds = sorted(self.inter_bonds + self.intra_bonds)

    def __repr__(self):
        return f"BipartiteLattice({self.nearest_bonds})"

    def _repr_html_(self):
        """
        HTML representation for visualizing the lattice with nodes and bonds.
        """
        fig, ax = plt.subplots(figsize=(3, 3))

        # Create coordinates for nodes
        coordinates = {
            idx: (x, self.ny - y - 1)
            for y in range(self.ny)
            for x, idx in enumerate(self.indices[y * self.nx: (y + 1) * self.nx])
        }

        # Draw bonds
        for i, j in self.nearest_bonds:
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[j]
            ax.plot([x1, x2], [y1, y2], ls="-", c="lightgray", lw=2)

        # Draw bonds
        for i, j in self.not_nearest_inter_bonds:
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[j]
            ax.plot([x1, x2], [y1, y2], ls="-", c="lightgray", lw=1)

        # Draw intra_bonds
        for i, j in self.intra_bonds:
            color = COLOR_A if i in self.indices_a else COLOR_B
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[j]
            ax.plot([x1, x2], [y1, y2], ls="-", c=color, lw=1)

        # Draw nodes
        for idx, (x, y) in coordinates.items():
            color = COLOR_A if idx in self.indices_a else COLOR_B
            ax.plot(x, y, 'o', color=color, markersize=10)
            ax.text(x, y, idx, c="w", ha="center", va="center", family="sans-serif", weight="medium")

        # Configure plot
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(-0.5, self.nx - 0.5)
        ax.set_ylim(-0.5, self.ny - 0.5)

        # Render plot as HTML
        buf = BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        data = base64.b64encode(buf.getvalue()).decode('utf-8')
        buf.close()
        return f"<img src='data:image/png;base64,{data}'/>"


