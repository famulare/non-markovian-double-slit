# Non-Markovian Double Slit

A pedagogical, computational double-slit model in Jacob Barandes's framework of
indivisible stochastic processes.

The implementation has a deliberate boundary:

1. Hamiltonian/Hilbert-space dynamics constructs a family of unistochastic laws,
   `Gamma[i, j] = abs(U[i, j])**2`.
2. The stochastic-process layer stores only nonnegative transition laws and
   probability distributions.
3. Event generation and the alternative path-space realizers use positive
   probabilities only.

This is not a claim that one-time probabilities obey a closed Markov evolution.
The failure of stochastic-matrix composition is the feature being demonstrated.

## Reproduce

Python environments and commands are managed with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run pytest
uv run python scripts/reproduce.py --profile quick
uv run python scripts/reproduce_leblond.py
```

The mathematical construction and interpretation are developed in
[`docs/derivation.md`](docs/derivation.md). The
[`prior-work audit`](docs/prior-work-audit.md) records exactly what is replicated,
corrected, and added so that original results are not conflated with prior work.

![Finite transverse law versus two-dimensional integrated current](docs/figures/exact_vs_2d.png)

## Scope

- Exact finite cyclic-grid construction of the tuple
  `(C, T, T0, Gamma, p, A)`.
- Direct versus artificially divided stochastic predictions.
- Both-slit, one-slit, phase-shift, event-buildup, partial which-path, coherent
  uncomputation, and conditional quantum-eraser demonstrations.
- Two inequivalent positive path-space realizers with the same observable
  marginals.
- A supporting two-dimensional finite-barrier simulation with a
  time-integrated probability-current screen.

## Repository map

- `src/non_markovian_double_slit/process.py`: positive-law tuple; no amplitude
  or Hamiltonian dependency.
- `src/non_markovian_double_slit/finite.py`: finite Hilbert-space construction
  of the law.
- `src/non_markovian_double_slit/realizers.py`: two positive path measures.
- `src/non_markovian_double_slit/which_path.py`: record overlap and physical
  erasure cases.
- `src/non_markovian_double_slit/lattice2d.py`: supporting split-step model.
- `scripts/`: reproducible figures and machine-readable diagnostics.
- `tests/`: semantic and numerical invariants.
- `literature/`: licensed source PDFs and provenance.

## License

Code is MIT licensed. Original prose and figures are CC BY 4.0. Third-party
papers in `literature/` retain their own licenses; see `THIRD_PARTY.md`.
