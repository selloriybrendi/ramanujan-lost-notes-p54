#!/usr/bin/env python3
"""audit_v3.py — mustaqil audit of Section 4 of the v2 (journal) / v3 (Zenodo) paper.

Siegel functions are implemented LITERALLY from the q-product (eq. (8) of the paper):
    g_a(tau) = -q^{B2(a1)/2} e^{pi i a2 (a1-1)} (1-q_z) prod_{n>=1} (1-q^n q_z)(1-q^n/q_z),   z = a1 tau + a2,
independently of the theta-function code.  Every intermediate claim of Lemmas 6-7, Theorem 8 and
Remarks 1-2 is then compared numerically (mpmath, 40 digits).  Prints MOS/FARQ per claim; exit 0 iff all MOS.
"""
import sys
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 40
OK = True


def rep(msg, ok):
    global OK
    OK = OK and bool(ok)
    print(("MOS  " if ok else "FARQ ") + msg)


def B2(x):
    return x * x - x + mp.mpf(1) / 6


def g_siegel(a1, a2, tau, N=4000):
    """Literal Siegel q-product, eq. (8). a1, a2 rational (mpf), tau in H."""
    q = mp.exp(2j * mp.pi * tau)
    z = a1 * tau + a2
    qz = mp.exp(2j * mp.pi * z)
    val = -mp.exp(2j * mp.pi * tau * B2(a1) / 2) * mp.exp(1j * mp.pi * a2 * (a1 - 1)) * (1 - qz)
    for n in range(1, N + 1):
        t = q ** n
        f1 = 1 - t * qz
        f2 = 1 - t / qz
        val *= f1 * f2
        if abs(t) < mp.mpf(10) ** -50:
            break
    return val


def prodq(q, res, N=15):
    p = mp.mpc(1)
    for r in res:
        k = 0
        while True:
            t = q ** (r + N * k)
            if abs(t) < mp.mpf(10) ** -50 or k > 100000:
                break
            p *= (1 - t)
            k += 1
    return p


def f(q, a, b):  # Ramanujan f(-q^a,-q^b) via Jacobi triple product
    return prodq(q, [a, b, a + b], N=a + b)


def u1(tau):
    q = mp.exp(2j * mp.pi * tau)
    return q * f(q, 2, 13) / f(q, 7, 8)


def u2(tau):
    q = mp.exp(2j * mp.pi * tau)
    return q * f(q, 1, 14) / f(q, 4, 11)


def act(M, tau):
    a, b, c, d = M
    return (a * tau + b) / (c * tau + d)


F15 = lambda r: mp.mpf(r) / 15
tau = mp.mpc('0.13', '0.21')
tau2 = mp.mpc('-0.31', '0.37')
tol = mp.mpf(10) ** -25

print("== Lemma 6: e(r) = 15 B2(r/15)/2 = r^2/30 - r/2 + 5/4 (aniq kasr)")
for r, want in [(2, Fr(23, 60)), (7, Fr(-37, 60)), (1, Fr(47, 60)), (4, Fr(-13, 60))]:
    e = Fr(15, 2) * (Fr(r, 15) ** 2 - Fr(r, 15) + Fr(1, 6))
    e2 = Fr(r * r, 30) - Fr(r, 2) + Fr(5, 4)
    rep(f"e({r}) = {e} = {e2} = {want}", e == want and e2 == want)

print("== Lemma 6: u1 = g_(2/15,0)(15 tau)/g_(7/15,0)(15 tau), u2 = g_(1/15,0)/g_(4/15,0)  (literal Siegel vs theta code)")
for t in (tau, tau2):
    tp = 15 * t
    lhs1 = g_siegel(F15(2), 0, tp) / g_siegel(F15(7), 0, tp)
    lhs2 = g_siegel(F15(1), 0, tp) / g_siegel(F15(4), 0, tp)
    rep(f"u1 at tau={mp.nstr(t, 4)}: |g2/g7 - u1| = {mp.nstr(abs(lhs1 - u1(t)), 3)}", abs(lhs1 - u1(t)) < tol)
    rep(f"u2 at tau={mp.nstr(t, 4)}: |g1/g4 - u2| = {mp.nstr(abs(lhs2 - u2(t)), 3)}", abs(lhs2 - u2(t)) < tol)

print("== Lemma 7 proof, step 1 (K1 for quotients, row-vector action): g_(1/15,0)(g' t')/g_(4/15,0)(g' t') = g_(7/15,6)(t')/g_(28/15,24)(t')")
Gp = (7, 90, 1, 13)
assert Gp[0] * Gp[3] - Gp[1] * Gp[2] == 1
for t in (tau, tau2):
    tp = 15 * t
    gtp = act(Gp, tp)
    lhs = g_siegel(F15(1), 0, gtp) / g_siegel(F15(4), 0, gtp)
    rhs = g_siegel(F15(7), 6, tp) / g_siegel(F15(28), 24, tp)
    rep(f"K1 quotient at tau={mp.nstr(t, 4)}: |LHS-RHS| = {mp.nstr(abs(lhs - rhs), 3)}", abs(lhs - rhs) < tol)
    rep("  15*gamma(tau) == gamma'(15 tau)", abs(15 * act((7, 6, 15, 13), t) - gtp) < tol)

print("== Lemma 7 proof, step 2 (K2 constants from the literal product)")
for t in (tau, tau2):
    tp = 15 * t
    c1 = g_siegel(F15(7), 6, tp) / g_siegel(F15(7), 0, tp)
    c2 = g_siegel(F15(28), 24, tp) / g_siegel(-F15(2), 0, tp)
    c3 = g_siegel(-F15(2), 0, tp) / g_siegel(F15(2), 0, tp)
    rep(f"g_(7/15,6)/g_(7/15,0) = e^(4 pi i/5): |diff| = {mp.nstr(abs(c1 - mp.exp(4j * mp.pi / 5)), 3)}", abs(c1 - mp.exp(4j * mp.pi / 5)) < tol)
    rep(f"g_(28/15,24)/g_(-2/15,0) = e^(-6 pi i/5): |diff| = {mp.nstr(abs(c2 - mp.exp(-6j * mp.pi / 5)), 3)}", abs(c2 - mp.exp(-6j * mp.pi / 5)) < tol)
    rep(f"g_(-2/15,0) = -g_(2/15,0): |diff| = {mp.nstr(abs(c3 + 1), 3)}", abs(c3 + 1) < tol)
rep("(-1)^(b1 b2 + b1 + b2) for b=(0,6) is +1 and for b=(2,24) is (-1)^74 = +1", (-1) ** (0 * 6 + 0 + 6) == 1 and (-1) ** (2 * 24 + 2 + 24) == 1 and 2 * 24 + 2 + 24 == 74)
rep("exponents: 6*7/15 = 14/5 -> e^(14 pi i/5) = e^(4 pi i/5);  24*2/15 = 16/5 -> e^(-16 pi i/5) = e^(-6 pi i/5)",
    Fr(6 * 7, 15) == Fr(14, 5) and Fr(14, 5) - 2 == Fr(4, 5) and Fr(24 * 2, 15) == Fr(16, 5) and Fr(-16, 5) + 2 == Fr(-6, 5))
rep("4/5 + 6/5 = 2  (so -e^(2 pi i) = -1)", Fr(4, 5) + Fr(6, 5) == 2)

print("== Lemma 7 statement: u2(gamma tau) u1(tau) = -1")
for t in (tau, tau2, mp.mpc('0.41', '0.07')):
    v = u2(act((7, 6, 15, 13), t)) * u1(t)
    rep(f"  tau={mp.nstr(t, 4)}: |u2(g t) u1(t) + 1| = {mp.nstr(abs(v + 1), 3)}", abs(v + 1) < tol)

print("== Theorem 8(i): rational bounds and the printed digits")
def u2_partial(q, K):
    p = q
    for k in range(K + 1):
        p *= (1 - q ** (15 * k + 1)) * (1 - q ** (15 * k + 14))
        p /= (1 - q ** (15 * k + 4)) * (1 - q ** (15 * k + 11))
    return p
P4 = u2_partial(Fr(-1, 2), 4)
P6 = u2_partial(Fr(-7, 10), 6)
S4 = 8 * Fr(1, 2) ** 76 / (1 - Fr(1, 2) ** 15)
S6 = 8 * Fr(7, 10) ** 106 / (1 - Fr(7, 10) ** 15)
rep(f"P_4(-1/2) = {float(P4):.17f} begins with -0.79954704982494", str(mp.mpf(P4.numerator) / P4.denominator).startswith("-0.79954704982494"))
rep(f"P_6(-7/10) = {float(P6):.17f} begins with -1.51863908249240", str(mp.mpf(P6.numerator) / P6.denominator).startswith("-1.51863908249240"))
rep(f"S(-1/2,K=4) = {float(S4):.3e} < 1.1e-22", S4 < Fr(11, 10) * Fr(1, 10 ** 22))
rep(f"S(-7/10,K=6) = {float(S6):.3e} < 3.1e-16", S6 < Fr(31, 10) * Fr(1, 10 ** 16))
rep("|q|^(15(K+1)+1) <= 1/2 for both (hypothesis of the tail bound)", Fr(1, 2) ** 76 <= Fr(1, 2) and Fr(7, 10) ** 106 <= Fr(1, 2))
rep("u2(-1/2) > -1 > u2(-7/10) (with tails)", P4 * (1 + S4 + S4 * S4) > -1 and P6 * (1 - S6) < -1)
rep("-log(1-x) <= 2x for x<=1/2 (used in the tail bound): check at x=1/2", -mp.log(mp.mpf(1) / 2) <= 1)

print("== Theorem 8: numerical constants printed in the paper")
y0 = mp.findroot(lambda y: (u2(mp.mpc(0.5, y)) + 1).real, (mp.mpf('0.08'), mp.mpf('0.1')), solver='bisect', tol=mp.mpf(10) ** -36)
tau0 = mp.mpc(0.5, y0)
q0 = mp.exp(2j * mp.pi * tau0).real
rep(f"q0 = {mp.nstr(q0, 20)} begins with -0.570655687750373", mp.nstr(q0, 16).startswith("-0.570655687750373"))
rep(f"y0 = {mp.nstr(y0, 20)} begins with 0.089281029042", mp.nstr(y0, 14).startswith("0.089281029042"))
rep("y0 = -log|q0|/(2 pi)", abs(y0 + mp.log(abs(q0)) / (2 * mp.pi)) < tol)
rep("-7/10 < q0 < -1/2", -0.7 < q0 < -0.5)
tau1_exact = act((13, -6, -15, 7), tau0)
rep("gamma^-1 = (13 -6; -15 7): gamma(gamma^-1 tau0) = tau0 exactly", abs(act((7, 6, 15, 13), tau1_exact) - tau0) < tol)
rep(f"Re(gamma^-1 tau0) = {mp.nstr(tau1_exact.real, 15)} == 0.117021434026 - 1 (the paper's 'mod 1')", abs(tau1_exact.real + 1 - mp.mpf('0.117021434026')) < 5e-13)
tau1 = mp.mpc(tau1_exact.real + 1, tau1_exact.imag)
q1 = mp.exp(2j * mp.pi * tau1)
rep("q(tau1) = q(gamma^-1 tau0): shift by 1 does not change q", abs(q1 - mp.exp(2j * mp.pi * tau1_exact)) < tol)
rep(f"tau1 = {mp.nstr(tau1.real, 15)} + {mp.nstr(tau1.imag, 15)} i rounds to 0.117021434026 + 0.043690294673 i",
    abs(tau1.real - mp.mpf('0.117021434026')) < 5e-13 and abs(tau1.imag - mp.mpf('0.043690294673')) < 5e-13)
rep(f"q(tau1) = {mp.nstr(q1, 15)} rounds to 0.563611486667 + 0.509757510046 i",
    abs(q1.real - mp.mpf('0.563611486667')) < 5e-13 and abs(q1.imag - mp.mpf('0.509757510046')) < 5e-13)
A1, B1 = f(q1, 7, 8), q1 * f(q1, 2, 13)
rep(f"A1(tau1) = B1(tau1): |A-B| = {mp.nstr(abs(A1 - B1), 3)}", abs(A1 - B1) < mp.mpf(10) ** -30)
rep("A1^3 - B1^3 = 0 and (A1-B1)^3 = 0 at tau1", abs(A1 ** 3 - B1 ** 3) < mp.mpf(10) ** -30)
A2v, B2v = f(mp.mpf(q0), 4, 11), mp.mpf(q0) * f(mp.mpf(q0), 1, 14)
rep(f"A2(q0) + B2(q0) = 0: |A+B| = {mp.nstr(abs(A2v + B2v), 3)}; A2^3+B2^3 = {mp.nstr(abs(A2v**3 + B2v**3), 3)}", abs(A2v + B2v) < mp.mpf(10) ** -30 and abs(A2v ** 3 + B2v ** 3) < mp.mpf(10) ** -30)

print("== Remark 1: the four one-term products are nonzero at tau0 and tau1 (and equal the cubes)")
rec = [("(8.2.6)", f(q1, 2, 3) * f(q1, 5, 10) / f(q1, 1, 4), A1 + B1), ("(8.2.8)", f(q1, 6, 9) * f(q1, 5, 10) ** 3 / f(q1, 3, 12), A1 ** 3 + B1 ** 3),
       ("(8.2.10)", f(mp.mpf(q0), 1, 4) * f(mp.mpf(q0), 5, 10) / f(mp.mpf(q0), 2, 3), A2v - B2v), ("(8.2.12)", f(mp.mpf(q0), 3, 12) * f(mp.mpf(q0), 5, 10) ** 3 / f(mp.mpf(q0), 6, 9), A2v ** 3 - B2v ** 3)]
for name, prod, direct in rec:
    rep(f"{name}: |product| = {mp.nstr(abs(prod), 6)} > 0 and equals the direct combination (|diff| = {mp.nstr(abs(prod - direct), 3)})", abs(prod) > 1e-6 and abs(prod - direct) < mp.mpf(10) ** -28)

print("== Remark 2: Gamma_0(15) elements (checked to 30 digits in the paper)")
def const(fn, M, pts=(tau, tau2)):
    vals = [fn(act(M, t), t) for t in pts]
    return vals[0] if abs(vals[0] - vals[1]) < mp.mpf(10) ** -25 else None
for M in [(1, 1, 15, 16), (16, 1, 15, 1)]:
    c = const(lambda gt, t: u1(gt) / u1(t), M)
    c2 = const(lambda gt, t: u2(gt) / u2(t), M)
    rep(f"{M}: u1(g t) = u1(t) and u2(g t) = u2(t)", c is not None and abs(c - 1) < tol and c2 is not None and abs(c2 - 1) < tol)
c = const(lambda gt, t: u1(gt) * u1(t), (4, 1, 15, 4))
rep("(4 1; 15 4): u1(g t) = 1/u1(t)", c is not None and abs(c - 1) < tol)
c = const(lambda gt, t: u1(gt) * u1(t), (11, 8, 15, 11))
rep("(11 8; 15 11): u1(g t) = 1/u1(t)  (upper-left == -4 mod 15)", c is not None and abs(c - 1) < tol)
c = const(lambda gt, t: u2(gt) / u1(t), (2, 1, 15, 8))
rep("(2 1; 15 8): u2(g t) = -u1(t)", c is not None and abs(c + 1) < tol)
c = const(lambda gt, t: u2(gt) / u1(t), (13, 6, 15, 7))
rep("(13 6; 15 7): u2(g t) = -u1(t)", c is not None and abs(c + 1) < tol)
c = const(lambda gt, t: u2(gt) * u1(t), (8, 1, 15, 2))
rep("(8 1; 15 2): u2(g t) = -1/u1(t)", c is not None and abs(c + 1) < tol)

print("== Section 5.1 remark: u1 = omega at tau ~ 0.29434 + 0.06065 i; cofactor relations")
w = mp.exp(2j * mp.pi / 3)
z = mp.findroot(lambda t: u1(t) - w, mp.mpc('0.2943', '0.0607'), tol=mp.mpf(10) ** -30)
rep(f"u1(tau) = omega at tau = {mp.nstr(z, 10)} (|u1-omega| = {mp.nstr(abs(u1(z) - w), 3)}); matches 0.29434+0.06065i", abs(u1(z) - w) < tol and abs(z - mp.mpc('0.29434', '0.06065')) < 1e-4)
rep("1 - u + u^2 = 0 iff u in {-omega, -omega^2}", abs(1 + w + w * w) < tol and abs(1 - (-w) + (-w) ** 2) < tol)

print("\nAUDIT:", "HAMMASI MOS" if OK else "FARQ BOR")
sys.exit(0 if OK else 1)
