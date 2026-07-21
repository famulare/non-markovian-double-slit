# Source notes: direct double-slit predecessors

## Barandes, Lecture Notes #22

[*Double-Slit Interference for a Unistochastic
Process*](https://shared.jacobbarandes.com/documents/double-slit-interference-unistochastic-lecture-notes)
is the direct two-configuration predecessor. It derives the coherent cross term,
shows loss of interference with a persistent orthogonal record, and discusses a
record that later "forgets."

The reviewed PDF states that it is all rights reserved, so it is linked rather
than copied into `literature/`. The local-forgetting operation needs the physical
qualification developed in `docs/prior-work-audit.md`.

## LeBlond, arXiv:2511.02863v1

David LeBlond extends the lecture-note argument to a physical electron example
with 2,000 sampled positions and includes its R implementation. Parameters,
kernel equations, and listings are reproduced in
`scripts/reproduce_leblond.py`. The source is CC0 and the PDF is archived in
`literature/`.

Key locators in the v1 PDF:

- assumptions and geometry: pp. 3-4;
- amplitudes and free kernel: Eqs. 1-11, pp. 6-8;
- results and stated fringe scale: pp. 8-9;
- R implementation: Listings 1-2, pp. 10-14.

Source: [arXiv:2511.02863](https://arxiv.org/abs/2511.02863).

## Bauer, arXiv:1304.1767v3

M. Bauer's spatial point-slit initial condition and free evolution appear in
Eqs. 5-9. The solution gives the relative phase and familiar fringe positions
analytically. The delta-function slit state is schematic and not
normalized as a finite vector, so this project uses finite-width normalized box
modes instead.

Source: [arXiv:1304.1767](https://arxiv.org/abs/1304.1767).
