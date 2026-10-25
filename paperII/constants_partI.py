"""Interval arithmetic for the constants of Part I of "The effective mass of the Pauli-Fierz electron" (S. Futakuchi, 2026).

Constants (sigma = K/2m):
  x_K/alpha = 4/(3 pi) log(1+sigma)
  t_K/alpha = 3 pi L(sigma)^2 / (16 log(1+sigma)),  L(s) = sqrt(s(1+s)) - arsinh(sqrt s)
  q_K       = 2 (log(1+sigma) - sigma/(1+sigma)) / log(1+sigma)
  rho_K     = 2 sqrt(q_K t_K) + 2 t_K = c1 sqrt(alpha) + c2 alpha
Also checks the two closed-form integrals of the appendix by quadrature:
  int_0^K k^2 g beta dk = 2 m^2 L(sigma)        (g = 1/Omega_k, beta^2 = (2m+k)/(4k))
  int_0^K 4 m^2 k/(2m+k)^2 dk = 4 m^2 (log(1+sigma) - sigma/(1+sigma))
The interval values are the constants used in the proof of the explicit corollary
of Part I (appendix on the constants); the quadratures are cross-checks only.
"""
import hashlib, io, json, platform, sys
import mpmath
from mpmath import iv, mp

iv.dps = 40
mp.dps = 40


def consts(sig):
    s = iv.mpf(sig)
    lg = iv.log(1 + s)
    L = iv.sqrt(s * (1 + s)) - iv.log(iv.sqrt(s) + iv.sqrt(1 + s))  # arsinh(sqrt s)
    x = 4 / (3 * iv.pi) * lg
    t = 3 * iv.pi * L ** 2 / (16 * lg)
    q = 2 * (lg - s / (1 + s)) / lg
    c1 = 2 * iv.sqrt(q * t)
    c2 = 2 * t
    return dict(x=x, two_x=2 * x, t=t, q=q, c1=c1, c2=c2, L=L)


def quad_checks(sig, m=1.0):
    K = 2 * m * sig
    g = lambda k: 1 / (k + k * k / (2 * m))
    beta = lambda k: mp.sqrt((2 * m + k) / (4 * k))
    I1 = mp.quad(lambda k: k * k * g(k) * beta(k), [0, K])
    L = mp.sqrt(sig * (1 + sig)) - mp.asinh(mp.sqrt(sig))
    I2 = mp.quad(lambda k: 4 * m * m * k / (2 * m + k) ** 2, [0, K])
    return dict(I1_minus_closed=float(abs(I1 - 2 * m * m * L)),
                I2_minus_closed=float(abs(I2 - 4 * m * m * (mp.log(1 + sig) - sig / (1 + sig)))))


def iv2(a):
    return [float(a.a), float(a.b)]


res = {}
c = consts(iv.mpf(1) / 2)
res["sigma_half"] = {k: iv2(v) for k, v in c.items()}
res["sigma_half_str"] = {k: str(v) for k, v in c.items()}
# safe rounding checks for the corollary
res["check_2x_ge_0.3441"] = bool(c["two_x"].a >= iv.mpf("0.3441").b)
res["check_c1_le_0.2985"] = bool(c["c1"].b <= iv.mpf("0.2985").a)
res["check_c2_le_0.1252"] = bool(c["c2"].b <= iv.mpf("0.1252").a)
a = iv.mpf(1) / iv.mpf("137.036")
hb = c["two_x"] * a / (1 + c["c1"] * iv.sqrt(a) + c["c2"] * a)
res["h_minus_1_at_alpha_1_137"] = iv2(hb)
res["saturation"] = iv2(1 + c["two_x"] / c["c2"])
res["quadrature"] = {s: quad_checks(mp.mpf(s), m=mp.mpf(1)) for s in ("0.5", "0.1", "2")}
src = open(__file__, "rb").read()
res["input_hash"] = hashlib.sha256(src).hexdigest()[:16]
res["python"] = sys.version.split()[0]
res["mpmath"] = mpmath.__version__
res["platform"] = platform.platform()
ok = res["check_2x_ge_0.3441"] and res["check_c1_le_0.2985"] and res["check_c2_le_0.1252"] \
    and all(v["I1_minus_closed"] < 1e-20 and v["I2_minus_closed"] < 1e-20 for v in res["quadrature"].values())
res["verdict"] = "PASS" if ok else "COUNTEREXAMPLE"
out = __file__.replace(".py", "_result.json")
with io.open(out, "w", encoding="utf-8", newline="\n") as fh:
    json.dump(res, fh, indent=1)
print(json.dumps(res, indent=1))
