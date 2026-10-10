#!/usr/bin/env python3
"""nol_sertifikat.py — rigorous certificate for Section "Why a note was necessary".

Part A (exact, rational arithmetic): the function u2(q) = B/A of Entry 8.2.2,
    u2(q) = q (q,q^14;q^15)_inf / (q^4,q^11;q^15)_inf,
satisfies u2(-1/2) > -1 and u2(-7/10) < -1.  Since u2 is continuous on (-1,0)
and A != 0 there, the intermediate value theorem gives q0 in (-7/10,-1/2) with
B(q0) = -A(q0), hence A^3 + B^3 = 0 at q0.  All bounds are proved: the partial
product P_K is an exact rational number, and the tail is enclosed by
|ln(tail)| <= S with S an explicit rational upper bound.

Part B (numerical, 40 digits, mpmath): the conjugation lemma
    u2(gamma tau) * u1(tau) = -1,   gamma = (7 6; 15 13) in Gamma_0(15),
the location of the real zero tau0 = 1/2 + i y0, and the image tau1 = gamma^{-1} tau0
at which A = B for Entry 8.2.1 (so A^3 - B^3 = 0).  Part B is a check of
statements that are proved in the paper; Part A is itself the proof.

Exit status 0 iff every certificate passes.
"""
import sys
from fractions import Fraction as Fr

OK = True


def rep(msg, ok):
    global OK
    OK = OK and ok
    print(("MOS  " if ok else "FARQ ") + msg)


# ---------------------------------------------------------------- Part A
def u2_partial(q, K):
    """Exact rational: q * prod_{k=0}^{K} (1-q^{15k+1})(1-q^{15k+14}) / ((1-q^{15k+4})(1-q^{15k+11}))."""
    p = q
    for k in range(K + 1):
        p *= (1 - q ** (15 * k + 1)) * (1 - q ** (15 * k + 14))
        p /= (1 - q ** (15 * k + 4)) * (1 - q ** (15 * k + 11))
    return p


def tail_log_bound(x, K):
    """Rational S with |ln prod_{k>K} r_k| <= S for |q| = x <= 1/2^(1/15) ... (we only use x <= 7/10).
    For |z| <= x^m <= 1/2:  |ln(1 - z)| <= -ln(1-|z|) <= 2|z|.  Each r_k has four factors, so
    |ln r_k| <= 2 (x^{15k+1} + x^{15k+4} + x^{15k+11} + x^{15k+14}) <= 8 x^{15k+1}, and
    sum_{k>K} 8 x^{15k+1} = 8 x^{15(K+1)+1} / (1 - x^15)."""
    assert x ** (15 * (K + 1) + 1) <= Fr(1, 2)
    return 8 * x ** (15 * (K + 1) + 1) / (1 - x ** 15)


def exp_upper(s):
    """Rational upper bound for e^s, 0 <= s <= 1:  e^s <= 1 + s + s^2  (true for 0<=s<=1)."""
    assert 0 <= s <= 1
    return 1 + s + s * s


def exp_lower(s):
    """Rational lower bound for e^{-s}, 0 <= s <= 1:  e^{-s} >= 1 - s."""
    assert 0 <= s <= 1
    return 1 - s


print("== Part A: exact rational certificate for the real zero of A^3+B^3 (Entry 8.2.2)")
for q, K, side in [(Fr(-1, 2), 4, "above"), (Fr(-7, 10), 6, "below")]:
    P = u2_partial(q, K)
    S = tail_log_bound(abs(q), K)
    # P < 0, so u2 = P e^{theta} with |theta| <= S lies in [P e^{S}, P e^{-S}]
    lo = P * exp_upper(S)
    hi = P * exp_lower(S)
    assert lo <= hi and P < 0
    print(f"   q = {q}: P_{K} = {float(P):.15f} (exact rational, denominator {len(str(P.denominator))} digits)")
    print(f"             tail bound S = {float(S):.3e};  u2(q) in [{float(lo):.15f}, {float(hi):.15f}]")
    if side == "above":
        rep(f"u2({q}) > -1  (lower bound {float(lo):.12f} > -1)", lo > -1)
    else:
        rep(f"u2({q}) < -1  (upper bound {float(hi):.12f} < -1)", hi < -1)
print("   => by the intermediate value theorem there is q0 in (-7/10, -1/2) with u2(q0) = -1,")
print("      i.e. B(q0) = -A(q0) and A^3 + B^3 = 0 at q0.  (A(q0) != 0: product of nonzero factors.)")

# ---------------------------------------------------------------- Part B
try:
    import mpmath as mp
except ImportError:
    print("mpmath yo'q — Part B o'tkazib yuborildi (Part A isbot uchun yetarli).")
    sys.exit(0 if OK else 1)

mp.mp.dps = 40


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


def f(q, a, b):  # Ramanujan f(-q^a,-q^b) by the Jacobi triple product
    return prodq(q, [a, b, a + b], N=a + b)


def u1(tau):
    q = mp.exp(2j * mp.pi * tau)
    return q * prodq(q, [2, 13]) / prodq(q, [7, 8])


def u2(tau):
    q = mp.exp(2j * mp.pi * tau)
    return q * prodq(q, [1, 14]) / prodq(q, [4, 11])


def act(M, tau):
    a, b, c, d = M
    return (a * tau + b) / (c * tau + d)


G = (7, 6, 15, 13)
assert G[0] * G[3] - G[1] * G[2] == 1

print("\n== Part B: conjugation lemma u2(gamma tau) u1(tau) = -1, gamma = (7 6; 15 13)")
worst = mp.mpf(0)
for tau in [mp.mpc('0.1', '0.3'), mp.mpc('0.37', '0.5'), mp.mpc('-0.2', '0.8'), mp.mpc('0.05', '0.12'),
            mp.mpc('0.6', '0.25'), mp.mpc('0.41', '0.07')]:
    e = u2(act(G, tau)) * u1(tau)
    worst = max(worst, abs(e + 1))
rep(f"max |u2(gamma tau) u1(tau) + 1| over 6 points = {mp.nstr(worst, 3)} < 1e-30", worst < mp.mpf(10) ** -30)

print("\n== Part B: the real zero tau0 = 1/2 + i y0 of Entry 8.2.2 and its image tau1 = gamma^{-1} tau0")
y0 = mp.findroot(lambda y: (u2(mp.mpc(0.5, y)) + 1).real, (mp.mpf('0.08'), mp.mpf('0.1')), solver='bisect',
                 tol=mp.mpf(10) ** -36)
tau0 = mp.mpc(0.5, y0)
q0 = mp.exp(2j * mp.pi * tau0)
print(f"   y0 = {mp.nstr(y0, 30)}")
print(f"   q0 = {mp.nstr(q0.real, 30)}   (Im q0 = {mp.nstr(q0.imag, 3)})")
A2, B2 = f(q0, 4, 11), q0 * f(q0, 1, 14)
rep(f"Entry 8.2.2 at q0: |A^3 + B^3| = {mp.nstr(abs(A2**3 + B2**3), 3)}", abs(A2 ** 3 + B2 ** 3) < mp.mpf(10) ** -30)
rep("q0 in (-7/10, -1/2)", -0.7 < q0.real < -0.5)

tau1 = act((13, -6, -15, 7), tau0)
tau1 = mp.mpc(tau1.real + 1, tau1.imag)  # shift into 0 <= Re < 1
q1 = mp.exp(2j * mp.pi * tau1)
A1, B1 = f(q1, 7, 8), q1 * f(q1, 2, 13)
print(f"   tau1 = {mp.nstr(tau1.real, 25)} + {mp.nstr(tau1.imag, 25)} i")
print(f"   q1   = {mp.nstr(q1, 25)}")
rep(f"Entry 8.2.1 at tau1: |A - B| = {mp.nstr(abs(A1 - B1), 3)}", abs(A1 - B1) < mp.mpf(10) ** -30)
rep(f"Entry 8.2.1 at tau1: |A^3 - B^3| = {mp.nstr(abs(A1**3 - B1**3), 3)}", abs(A1 ** 3 - B1 ** 3) < mp.mpf(10) ** -30)
rep(f"|A(tau1)| = {mp.nstr(abs(A1), 6)} != 0", abs(A1) > mp.mpf('0.01'))

print("\n== Part B: control — the recorded combinations are products, hence nonzero at tau0, tau1")
rec = {
    "(8.2.8)  A^3+B^3 (8.2.1) at tau1": f(q1, 6, 9) * f(q1, 5, 10) ** 3 / f(q1, 3, 12),
    "(8.2.12) A^3-B^3 (8.2.2) at tau0": f(q0, 3, 12) * f(q0, 5, 10) ** 3 / f(q0, 6, 9),
    "(8.2.6)  A+B (8.2.1) at tau1": f(q1, 2, 3) * f(q1, 5, 10) / f(q1, 1, 4),
    "(8.2.10) A-B (8.2.2) at tau0": f(q0, 1, 4) * f(q0, 5, 10) / f(q0, 2, 3),
}
for name, v in rec.items():
    rep(f"{name}: |value| = {mp.nstr(abs(v), 6)} > 0", abs(v) > mp.mpf('1e-6'))
# cross-check (8.2.8) and (8.2.12) against the direct cubes at the same points
rep("(8.2.8) matches A^3+B^3 at tau1", abs(rec["(8.2.8)  A^3+B^3 (8.2.1) at tau1"] - (A1 ** 3 + B1 ** 3)) < mp.mpf(10) ** -28)
rep("(8.2.12) matches A^3-B^3 at tau0", abs(rec["(8.2.12) A^3-B^3 (8.2.2) at tau0"] - (A2 ** 3 - B2 ** 3)) < mp.mpf(10) ** -28)

print("\nNATIJA:", "HAMMASI MOS" if OK else "FARQ BOR")
sys.exit(0 if OK else 1)
