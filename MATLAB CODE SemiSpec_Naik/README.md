# SemiSpec Toolkit

MATLAB / GNU Octave toolkit accompanying the paper:

> **SemiSpec: Eigenvalue- and Laplacian-Based Topological Indices of Semigraphs, Resistance-Distance Formulas and the Robustness of Semigraph Network Models**
> Sarojini Naik, Meenal M. Kaliwal and Shankar Naik

The toolkit computes eleven spectral and resistance-type invariants (spectral radius, energy, Estrada index, natural connectivity, Laplacian / signless-Laplacian / Laplacian-Estrada energies, algebraic connectivity, Laplacian spectral radius, Kirchhoff index, number of spanning trees) of the adjacency graph G_a(S) and the consecutive-adjacency graph G_c(S) of a semigraph S. It verifies every closed form in the paper to machine precision.

## Requirements

- MATLAB R2016b or later, **or** GNU Octave 6 or later
- No toolboxes required

## Quick start

```matlab
>> SemiSpec_Toolkit
```

Running the main file:
- verifies each closed form of Sections 4–6 of the paper
- writes Tables 1–4 (and Table 6, correlations) as CSV files
- redraws Figures 1–7
- generates the 400-semigraph random ensemble of Section 7 (`rng(2026)`)

## Entering your own semigraph

A semigraph is entered as a cell array of ordered s-edges:

```matlab
edges = {[1 2 3 4], [4 5 6], [2 6]};
check_semigraph(edges);          % enforces the semigraph axioms
[Aa, Ac] = sg_build(edges);      % adjacency and consecutive-adjacency matrices
idx = sg_indices(edges);         % all seventeen indices of Table 1
```

## Repository contents

| File / folder | Description |
|---|---|
| `SemiSpec_Toolkit.m` | Main MATLAB/Octave toolkit |
| `results/` | CSV tables produced by the toolkit |
| `figures/` | Figures 1–7 produced by the toolkit |

## How to cite

If you use this code, please cite the paper:

```
S. Naik, M. M. Kaliwal, S. Naik, SemiSpec: Eigenvalue- and Laplacian-Based Topological
Indices of Semigraphs, Resistance-Distance Formulas and the Robustness of Semigraph
Network Models, [Journal], [Year].
```

## Contact

Shankar Naik (corresponding author), Department of Mathematics, Anuvartik Mirji Bharatesh Institute of Technology (AMBIT), Belagavi, Karnataka, India — cme.shank@gmail.com

## License

Released under the MIT License (see `LICENSE`).
