# paperII

Ancillary code for paper [II], *The effective mass of the Pauli–Fierz electron: a cutoff-independent lower bound and the breakdown of mass renormalization*, [arXiv:2610.34125](https://arxiv.org/abs/2610.34125).

- `constants_partI.py`: interval-arithmetic evaluation (`mpmath.iv`, 40 digits) of the constants of Part I of [II] (appendix "Numerical values of the constants"): $x_K/\alpha$, $t_K/\alpha$, $q_K$ and $\rho_K=c_1\sqrt\alpha+c_2\alpha$ for $K=m$, and checks of the two closed-form integrals by quadrature. It checks closed-form arithmetic only and is not used in any proof.
- `constants_partI_result.json`: the output (verdict `PASS`).

Run with `python3 constants_partI.py` (tested with Python 3.9.6 and mpmath 1.3.0); the script writes `constants_partI_result.json` next to itself.
