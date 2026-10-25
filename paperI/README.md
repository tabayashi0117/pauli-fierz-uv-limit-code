# paperI

Ancillary code for paper [I], *The sharp leading asymptotics of the self-energy of the spinless Pauli–Fierz model*.

- `constants_interval.py`: rigorous enclosures, in interval arithmetic (`mpmath.iv`, 200 bits), of the constants of Appendix B of [I]: the first positive root $x_1$ of $\tan x=x$, $F[1]$ (by two closed forms, whose enclosures are checked to overlap), $c_e$ and $c_\alpha$, the crossover constants $2\pi c_\alpha$ and $8\pi^2c_e$, and (as a plain floating-point value) the crossover cutoff $\Lambda_c$ at $\alpha=1/137.036$. The printed enclosures are rounded outward.

Run with `python3 constants_interval.py` (requires `mpmath`; tested with Python 3.9.6 and mpmath 1.3.0).
