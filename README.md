# pauli-fierz-uv-limit-code

Ancillary code for two companion papers on the spinless Pauli–Fierz model with a sharp ultraviolet cutoff $\Lambda$ (one non-relativistic electron coupled to the quantized transverse electromagnetic field, no external potential):

- **[I]** S. Futakuchi, *The sharp leading asymptotics of the self-energy of the spinless Pauli–Fierz model*, [arXiv:2610.34122](https://arxiv.org/abs/2610.34122) [math-ph] (2026).
- **[II]** S. Futakuchi, *The effective mass of the Pauli–Fierz electron: a cutoff-independent lower bound and the breakdown of mass renormalization*, [arXiv:2610.34125](https://arxiv.org/abs/2610.34125) [math-ph] (2026).

The same scripts are included as ancillary files (`anc/`) of the arXiv submissions.

## Contents

| Directory | Paper | Contents |
|---|---|---|
| `paperI/` | [I] | `constants_interval.py`: rigorous interval-arithmetic enclosures (`mpmath.iv`) of the explicit constants of Appendix B of [I], e.g. $F[1]$ and the coefficient $c_\alpha$ of $\alpha^{2/7}\Lambda^{12/7}$ |
| `paperII/` | [II] | `constants_partI.py`: interval-arithmetic evaluation of the constants of Part I of [II] and quadrature checks of two closed-form integrals; `constants_partI_result.json`: its output |

Each directory has a `README.md` with the command line and the tested library versions.

## Role of the code

The theorems of both papers are proved analytically. The scripts evaluate explicit closed-form constants (with outward-rounded interval arithmetic where stated) and perform cross-checks; no proof depends on a floating-point computation.

## Requirements

Python 3 and `mpmath` (tested with Python 3.9.6 and mpmath 1.3.0). Each script runs in less than a second.

## Use of AI tools

The papers and this code were prepared with substantial assistance from Claude models (Anthropic), accessed through the Claude Code agent; see the disclosure section of each paper. The author takes full responsibility for the content.

## License

MIT License; see `LICENSE`.
