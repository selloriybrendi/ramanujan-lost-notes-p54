# Author: Otakhon U. Kenjaev (ORCID 0009-0009-3566-9285) · doi:10.5281/zenodo.23030496 · CC BY 4.0
# verifier.py — yagona sertifikat-skript
# "The two lost notes on page 54 of Ramanujan's lost notebook" maqolasining
# BARCHA ayniyatlarini (Entry 8.2.1: 8.2.6-8.2.9, Entry 8.2.2: 8.2.10-8.2.13,
# maqola v2 dagi Theorem 1-2 (circular summation, Borwein a(q^5)) va Theorem 3-4 (elementar) 400 q-hadgacha qayta tekshiradi.
#
# Texnika: y = q^(1/3) o'zgaruvchisida ishlaymiz (K y-had = K/3 q-had),
# bo'lish YO'Q — har kasrli ayniyat ko'paytirib (cross-multiply) sinaladi.
# Barcha koeffitsiyentlar butun son: hech qanday yaqinlashish/float yo'q.
#
# Ishlatish: python3 verifier.py [Nq]   (standart Nq=400)
import sys
import time

Nq = int(sys.argv[1]) if len(sys.argv) > 1 else 400
K = 3 * Nq + 3  # y-hadlar soni

def theta(a, b):
    """f(-y^a, -y^b) = sum_{n} (-1)^n y^{a n(n+1)/2 + b n(n-1)/2} (Jacobi uchlik mahsulotsiz, to'g'ridan qator)."""
    s = [0] * K
    n = 0
    while True:
        did = False
        for m in ((n, -n) if n else (0,)):
            e = a * m * (m + 1) // 2 + b * m * (m - 1) // 2
            if 0 <= e < K:
                s[e] += (-1) ** (m % 2)
                did = True
        if not did and n > 0:
            break
        n += 1
    return s

def kop(a, b):
    r = [0] * K
    for i, x in enumerate(a):
        if x:
            lim = K - i
            for j, y in enumerate(b[:lim]):
                if y:
                    r[i + j] += x * y
    return r

def kub(a):
    return kop(kop(a, a), a)

def sil(s, k):
    """y^k ga ko'paytirish (k>=0)."""
    r = [0] * K
    for i, x in enumerate(s):
        if x and i + k < K:
            r[i + k] = x
    return r

def ayr(a, b): return [x - y for x, y in zip(a, b)]
def qosh(a, b): return [x + y for x, y in zip(a, b)]
def skal(c, a): return [c * x for x in a]

nat = []
def sina(nom, chap, ong):
    ok = chap == ong
    if not ok:
        d = next(i for i in range(K) if chap[i] != ong[i])
        print(f"  {nom}: FARQ ✗  (birinchi farq y^{d}: {chap[d]} vs {ong[d]})")
    else:
        print(f"  {nom}: MOS ✓  ({Nq} q-had)")
    nat.append(ok)

t0 = time.time()
# ---- funksiyalar (y-dunyoda; q = y^3) ----
f14  = theta(3, 12)    # f(-q, -q^4)
f23  = theta(6, 9)     # f(-q^2, -q^3)
f312 = theta(9, 36)    # f(-q^3, -q^12)
f69  = theta(18, 27)   # f(-q^6, -q^9)
E5   = theta(15, 30)   # f(-q^5) = f(-q^5, -q^10)
f13_43 = theta(1, 4)   # f(-q^(1/3), -q^(4/3))
f23_1  = theta(2, 3)   # f(-q^(2/3), -q)

# Entry 8.2.1: A1 = f(-q^7, -q^8), B1 = q f(-q^2, -q^13)
A1 = theta(21, 24)
B1 = sil(theta(6, 39), 3)
# Entry 8.2.2: A2 = f(-q^4, -q^11), B2 = q f(-q, -q^14)
A2 = theta(12, 33)
B2 = sil(theta(3, 42), 3)

AB1 = kop(A1, B1)
AB2 = kop(A2, B2)
S1 = qosh(A1, B1)
D1 = ayr(A1, B1)
S2 = qosh(A2, B2)
D2 = ayr(A2, B2)

print("== Entry 8.2.1: A = f(-q^7,-q^8), B = q f(-q^2,-q^13)")
# (8.2.6) A+B = f23 E5 / f14   ->  (A+B) f14 == f23 E5
sina("(8.2.6)  (A+B)·f(-q,-q^4) = f(-q^2,-q^3)·f(-q^5)", kop(S1, f14), kop(f23, E5))
# (8.2.7) A-B = f(-q^(2/3),-q) + q^(2/3) f(-q^3,-q^12)
sina("(8.2.7)  A-B = f(-q^(2/3),-q) + q^(2/3)·f(-q^3,-q^12)", D1, qosh(f23_1, sil(f312, 2)))
# (8.2.8) A^3+B^3 = f69 E5^3 / f312  ->  (A^3+B^3) f312 == f69 E5^3
sina("(8.2.8)  (A^3+B^3)·f(-q^3,-q^12) = f(-q^6,-q^9)·f^3(-q^5)",
     kop(ayr(kub(A1), skal(-1, kub(B1))), f312), kop(f69, kub(E5)))
# (8.2.9) (A-B)^3 = f23 f14^3 / f312 + q^2 f312^3  ->  (A-B)^3 f312 == f23 f14^3 + q^2 f312^4
sina("(8.2.9)  (A-B)^3·f(-q^3,-q^12) = f(-q^2,-q^3)·f^3(-q,-q^4) + q^2·f^4(-q^3,-q^12)",
     kop(kub(D1), f312), qosh(kop(f23, kub(f14)), sil(kop(f312, kub(f312)), 6)))
# TEOREMA 1: A^3 - B^3 = (A-B)^3 + 3AB(A-B), komponentlar 8.2.9 va 8.2.7 bilan:
#   (A^3-B^3) f312 == [f23 f14^3 + q^2 f312^4] + 3AB(A-B) f312
T1chap = kop(ayr(kub(A1), kub(B1)), f312)
T1ong = qosh(qosh(kop(f23, kub(f14)), sil(kop(f312, kub(f312)), 6)),
             kop(skal(3, AB1), kop(qosh(f23_1, sil(f312, 2)), f312)))
sina("THEOREM 3 (elementar)  A^3-B^3 = (8.2.9) + 3AB·(8.2.7)", T1chap, T1ong)

print("== Entry 8.2.2: A = f(-q^4,-q^11), B = q f(-q,-q^14)")
# (8.2.10) A-B = f14 E5 / f23  ->  (A-B) f23 == f14 E5
sina("(8.2.10) (A-B)·f(-q^2,-q^3) = f(-q,-q^4)·f(-q^5)", kop(D2, f23), kop(f14, E5))
# (8.2.11) A+B = -q^(-1/3){ f(-q^(1/3),-q^(4/3)) - f69 }  ->  y(A+B) == f69 - f13_43
sina("(8.2.11) q^(1/3)·(A+B) = f(-q^6,-q^9) - f(-q^(1/3),-q^(4/3))", sil(S2, 1), ayr(f69, f13_43))
# (8.2.12) A^3-B^3 = f312 E5^3 / f69  ->  (A^3-B^3) f69 == f312 E5^3
sina("(8.2.12) (A^3-B^3)·f(-q^6,-q^9) = f(-q^3,-q^12)·f^3(-q^5)",
     kop(ayr(kub(A2), kub(B2)), f69), kop(f312, kub(E5)))
# (8.2.13) (A+B)^3 = -(1/q){ f14 f23^3 / f69 - f69^3 }  ->  q(A+B)^3 f69 == f69^4 - f14 f23^3
sina("(8.2.13) q·(A+B)^3·f(-q^6,-q^9) = f^4(-q^6,-q^9) - f(-q,-q^4)·f^3(-q^2,-q^3)",
     kop(sil(kub(S2), 3), f69), ayr(kop(f69, kub(f69)), kop(f14, kub(f23))))
# TEOREMA 2: A^3 + B^3 = (A+B)^3 - 3AB(A+B), komponentlar 8.2.13 va 8.2.11 bilan:
#   q(A^3+B^3) f69 == [f69^4 - f14 f23^3] - 3AB·q(A+B) f69, bunda q(A+B) = y^2·(f69-f13_43)/y^3...
#   aniq: q(A+B) = y^2 · [y(A+B)] = y^2 (f69 - f13_43)
T2chap = kop(sil(qosh(kub(A2), kub(B2)), 3), f69)
T2ong = ayr(ayr(kop(f69, kub(f69)), kop(f14, kub(f23))),
            kop(skal(3, AB2), kop(sil(ayr(f69, f13_43), 2), f69)))
sina("THEOREM 4 (elementar)  A^3+B^3 = (8.2.13) - 3AB·(8.2.11)", T2chap, T2ong)

# ---- v2: Entry 8.2.3 (circular summation) + Borwein a(q) ----
def borwein_a(k):
    """a(y^k) = sum_{m,n} y^{k(m^2+mn+n^2)}; m^2+mn+n^2 >= (3/4)max(|m|,|n|)^2 => |m|,|n| <= sqrt(4K/(3k))."""
    s = [0] * K
    M = int((4 * K / (3 * k)) ** 0.5) + 3
    for m in range(-M, M + 1):
        for n in range(-M, M + 1):
            e = k * (m * m + m * n + n * n)
            if e < K:
                s[e] += 1
    return s

a5 = borwein_a(15)  # a(q^5), q = y^3
print("== v2: circular summation (Entry 8.2.3) va Borwein a(q^5)")
# TEOREMA 3: A1^3 - B1^3 = q^2 f^3(-q^3,-q^12) + f(-q^2,-q^3) a(q^5)
sina("THEOREM 1 (circular)  A^3-B^3 = q^2 f^3(-q^3,-q^12) + f(-q^2,-q^3)·a(q^5)",
     ayr(kub(A1), kub(B1)), qosh(sil(kub(f312), 6), kop(f23, a5)))
# TEOREMA 4: q(A2^3 + B2^3) = f^3(-q^6,-q^9) - f(-q,-q^4) a(q^5)
sina("THEOREM 2 (circular)  q(A^3+B^3) = f^3(-q^6,-q^9) - f(-q,-q^4)·a(q^5)",
     sil(qosh(kub(A2), kub(B2)), 3), ayr(kub(f69), kop(f14, a5)))

print(f"\nJAMI: {sum(nat)}/{len(nat)} MOS · {time.time()-t0:.1f}s · EXIT={'0' if all(nat) else '1'}")
sys.exit(0 if all(nat) else 1)
