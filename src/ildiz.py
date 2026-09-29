# ildiz.py — ILDIZ darajasidagi tekshiruv (hech bir bosqich boshqasiga tayanmaydi):
#  (R1) theta-qator generatori o'zi to'g'rimi: f(-q^a,-q^b) YIG'INDI == Jacobi uchlik KO'PAYTMA (q^a;q^N)(q^b;q^N)(q^N;q^N)
#  (R2) Entry 8.2.3 UMUMIY holda: a,b — mustaqil formal o'zgaruvchilar, butun sonli ikki o'zgaruvchili ko'phad sifatida
#        f^3(ab^2,a^2b) + a f^3(b,a^3b^2) + b f^3(a,a^2b^3) == f(a,b) * A(ab),  A(x)=sum x^{m^2+mn+n^2}
#  (R3) Borwein kubik ayniyati a^3 = b^3 + c^3, bo'lishsiz: a(q)^3 (q;q)^3 (q^3;q^3)^3 == (q;q)^12 + 27 q (q^3;q^3)^12
#  (R4) Theorem 1-2 ni R2 ning xususiy holi sifatida emas, to'g'ridan-to'g'ri (bo'lishsiz) — mustaqil kod bilan
import sys
import time

t0 = time.time()
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
D = int(sys.argv[2]) if len(sys.argv) > 2 else 90
nat = []


def say(nom, ok, izoh=""):
    nat.append(ok)
    print(f"  {'MOS ✓' if ok else 'FARQ ✗'}  {nom} {izoh}")


def mul(a, b, n=N):
    r = [0] * (n + 1)
    for i, x in enumerate(a[: n + 1]):
        if x:
            for j, y in enumerate(b[: n + 1 - i]):
                if y:
                    r[i + j] += x * y
    return r


def sumtheta(a, b, n=N):  # f(-q^a,-q^b) yig'indi ta'rifidan
    s = [0] * (n + 1)
    k = 0
    while True:
        did = False
        for m in (k, -k) if k else (0,):
            e = a * m * (m + 1) // 2 + b * m * (m - 1) // 2
            if 0 <= e <= n:
                s[e] += (-1) ** (m % 2)
                did = True
        if not did and k > 0:
            break
        k += 1
    return s


def poch(a, M, n=N):  # (q^a;q^M)_inf ko'paytmadan
    r = [0] * (n + 1)
    r[0] = 1
    k = a
    while k <= n:
        r = [r[i] - (r[i - k] if i >= k else 0) for i in range(n + 1)]
        k += M
    return r


print(f"== R1: theta generatori == Jacobi ko'paytma ({N} had)")
for a, b in [
    (1, 4),
    (2, 3),
    (3, 12),
    (6, 9),
    (7, 8),
    (2, 13),
    (4, 11),
    (1, 14),
    (5, 10),
    (1, 2),
    (3, 6),
    (15, 30),
]:
    Mm = a + b
    say(
        f"f(-q^{a},-q^{b})",
        sumtheta(a, b) == mul(mul(poch(a, Mm), poch(b, Mm)), poch(Mm, Mm)),
    )
print(f"== R2: Entry 8.2.3 umumiy (a,b mustaqil), barcha a^i b^j, i+j <= {D}")


def ftheta2(p, r, s, t):  # f(a^p b^r, a^s b^t) — ikki o'zgaruvchili, butun daraja >=0
    d = {}
    k = 0
    while True:
        did = False
        for m in (k, -k) if k else (0,):
            u, v = m * (m + 1) // 2, m * (m - 1) // 2
            i, j = p * u + s * v, r * u + t * v
            if i + j <= D:
                d[(i, j)] = d.get((i, j), 0) + 1
                did = True
        if not did and k > 0:
            break
        k += 1
    return d


def m2(x, y):
    r = {}
    for (i1, j1), c1 in x.items():
        for (i2, j2), c2 in y.items():
            if i1 + i2 + j1 + j2 <= D:
                r[(i1 + i2, j1 + j2)] = r.get((i1 + i2, j1 + j2), 0) + c1 * c2
    return r


def shift(x, di, dj):
    return {(i + di, j + dj): c for (i, j), c in x.items() if i + di + j + dj <= D}


def add(*xs):
    r = {}
    for x in xs:
        for k, c in x.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c}


def cube(x):
    return m2(m2(x, x), x)


L = add(
    cube(ftheta2(1, 2, 2, 1)),
    shift(cube(ftheta2(0, 1, 3, 2)), 1, 0),
    shift(cube(ftheta2(1, 0, 2, 3)), 0, 1),
)
Aab = {}
M = int((4 * D / 6) ** 0.5) + 3
for m in range(-M, M + 1):
    for n in range(-M, M + 1):
        e = m * m + m * n + n * n
        if 2 * e <= D:
            Aab[(e, e)] = Aab.get((e, e), 0) + 1
R = add(m2(ftheta2(1, 0, 0, 1), Aab))
say(
    "8.2.3 bivariat",
    L == R,
    f"({len(L)} monomial, eng katta koef {max(abs(c) for c in L.values())})",
)
Rb = dict(R)
Rb[(3, 3)] = Rb.get((3, 3), 0) + 1
say("NAZORAT (buzilgan R2) — FARQ kutiladi", L != add(Rb), "")
print(f"== R3: Borwein a^3 = b^3 + c^3 (bo'lishsiz, {N} had)")
aq = [0] * (N + 1)
M = int((4 * N / 3) ** 0.5) + 3
for m in range(-M, M + 1):
    for n in range(-M, M + 1):
        e = m * m + m * n + n * n
        if e <= N:
            aq[e] += 1
e1, e3 = poch(1, 1), poch(3, 3)
e1c, e3c = mul(mul(e1, e1), e1), mul(mul(e3, e3), e3)
Lb = mul(mul(mul(aq, aq), aq), mul(e1c, e3c))
Rb3 = [
    x + 27 * y
    for x, y in zip(
        mul(mul(e1c, e1c), mul(e1c, e1c)), [0] + mul(mul(e3c, e3c), mul(e3c, e3c))[:N]
    )
]
say("a^3 (q;q)^3 (q^3;q^3)^3 == (q;q)^12 + 27q (q^3;q^3)^12", Lb == Rb3)
print(
    f"== R4: Theorem 1-2 to'g'ridan (JTP-ko'paytmadan qurilgan theta bilan, {N // 5} had)"
)
n5 = N // 5


def P(a, b):
    return mul(
        mul(poch(a, a + b, n5), poch(b, a + b, n5), n5), poch(a + b, a + b, n5), n5
    )


A1 = P(7, 8)
B1 = [0] + P(2, 13)[:n5]
A2 = P(4, 11)
B2 = [0] + P(1, 14)[:n5]
a5 = [aq[k // 5] if k % 5 == 0 else 0 for k in range(n5 + 1)]


def c(x):
    return mul(mul(x, x, n5), x, n5)


T1L = [x - y for x, y in zip(c(A1), c(B1))]
T1R = [x + y for x, y in zip([0, 0] + c(P(3, 12))[: n5 - 1], mul(P(2, 3), a5, n5))]
say("Theorem 1", T1L == T1R)
T2L = [0] + [x + y for x, y in zip(c(A2), c(B2))][:n5]
T2R = [x - y for x, y in zip(c(P(6, 9)), mul(P(1, 4), a5, n5))]
say("Theorem 2", T2L == T2R)
print(
    f"\nILDIZ: {sum(nat)}/{len(nat)} MOS · {time.time() - t0:.1f}s · EXIT={'0' if all(nat) else '1'}"
)
