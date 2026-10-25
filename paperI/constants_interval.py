#!/usr/bin/env python3
"""Rigorous enclosures of the constants of the paper (interval arithmetic, mpmath.iv).

Computes
  x_1    : first positive root of tan x = x (in (pi, 3pi/2)),
  F[1]   : (7/3)(3/10)^{3/7}(3pi/10)^{2/7} x_1^{6/7}   (Proposition B.1(b), first form),
  F[1]'  : (7/4)(2/3)^{3/7} C_3^{-3/7},  C_3 = 5^{5/3} / (2^{2/3} 3 x_1^2 (4pi/3)^{2/3})  (second form),
  c_e    : (6 pi^2)^{-4/7} F[1]            (coefficient of e^{4/7} Lambda^{12/7}),
  c_alpha: (4 pi)^{2/7} c_e                (coefficient of alpha^{2/7} Lambda^{12/7}, alpha = e^2/(4 pi)),
  2 pi c_alpha = 6.2416... (crossover constant of Remark 3.6, in terms of alpha); 8 pi^2 c_e = 38.058... (the same crossover in terms of e),
and the physical crossover cutoff Lambda_c = (2 pi c_alpha)^{7/2} alpha^{-5/2} at alpha = 1/137.036
(a plain floating-point value, not an enclosure).

Usage: python3 constants_interval.py      (requires mpmath >= 1.0; runs in < 1 s)
"""
import mpmath
from mpmath import iv, mp

PREC_BITS = 200
iv.prec = PREC_BITS
mp.prec = PREC_BITS


def g(x):
    """g(x) = sin x - x cos x; tan x = x on (pi, 3pi/2) iff g(x) = 0."""
    return iv.sin(x) - x * iv.cos(x)


def root_tan_x_equals_x():
    # bisection with interval sign checks on [4.4, 4.6], a subinterval of (pi, 3pi/2):
    # g'(x) = x sin x, negative on (pi, 3pi/2), so g is decreasing there; g(pi) = pi > 0, g(3pi/2) = -1 < 0.
    lo = mp.mpf(4.4)
    hi = mp.mpf(4.6)
    assert g(iv.mpf(lo)).a > 0 and g(iv.mpf(hi)).b < 0
    for _ in range(PREC_BITS - 10):
        mid = (lo + hi) / 2
        gm = g(iv.mpf(mid))
        if gm.a > 0:
            lo = mid
        elif gm.b < 0:
            hi = mid
        else:
            break
    return iv.mpf([lo, hi])


def main():
    x1 = root_tan_x_equals_x()
    pi = iv.pi
    F1 = iv.mpf(7) / 3 * (iv.mpf(3) / 10) ** (iv.mpf(3) / 7) * (3 * pi / 10) ** (iv.mpf(2) / 7) * x1 ** (iv.mpf(6) / 7)
    C3 = iv.mpf(5) ** (iv.mpf(5) / 3) / (iv.mpf(2) ** (iv.mpf(2) / 3) * 3 * x1 ** 2 * (4 * pi / 3) ** (iv.mpf(2) / 3))
    F1b = iv.mpf(7) / 4 * (iv.mpf(2) / 3) ** (iv.mpf(3) / 7) * C3 ** (-iv.mpf(3) / 7)
    ce = (6 * pi ** 2) ** (-iv.mpf(4) / 7) * F1
    ca = (4 * pi) ** (iv.mpf(2) / 7) * ce
    cross = 8 * pi ** 2 * ce
    nd = 30  # digits after the decimal point; endpoints are rounded outward by encl()

    def encl(x):
        s = mp.mpf(10) ** nd
        lo = mp.floor(mp.mpf(x.a) * s) / s
        hi = mp.ceil(mp.mpf(x.b) * s) / s
        f = lambda v: mp.nstr(v, nd + 5, min_fixed=-mp.inf, max_fixed=mp.inf)
        return "[" + f(lo) + ", " + f(hi) + "]"

    print("x_1          in", encl(x1))
    print("F[1]         in", encl(F1))
    print("F[1] (C_3)   in", encl(F1b))
    print("c_e          in", encl(ce))
    print("c_alpha      in", encl(ca))
    print("8 pi^2 c_e   in", encl(cross), "   2 pi c_alpha in", encl(2 * pi * ca))
    a_fs = mp.mpf(1) / mp.mpf("137.036")
    lam_c = (2 * mp.pi * ca.mid) ** (mp.mpf(7) / 2) * a_fs ** (-mp.mpf(5) / 2)
    print("Lambda_c (alpha = 1/137.036, value) =", mpmath.nstr(lam_c, 6))
    ok = (F1.a <= F1b.b) and (F1b.a <= F1.b)
    print("two closed forms of F[1] overlap:", ok)


if __name__ == "__main__":
    main()
