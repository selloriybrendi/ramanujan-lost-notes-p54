# Yo'qolgan izohlarni tiklash: T = c1*q^a1*M1 + c2*q^a2*M2 shaklini qidirish
# Baza: f14=f(-q,-q^4), f23=f(-q^2,-q^3), f312=f(-q^3,-q^12), f69=f(-q^6,-q^9), E5=f(-q^5)
# Monomial: darajalar yig'indisi 3 (kvotiyentga ruxsat), c=±1, q-siljish 0..2
from fractions import Fraction
N = 260
def theta(c4, c11):
    s = [0]*(N+1); n = 0
    while True:
        did = False
        for m in (n, -n) if n else (0,):
            e = c4*m*(m+1)//2 + c11*m*(m-1)//2
            if 0 <= e <= N: s[e] += (-1)**(m % 2); did = True
        if not did and n > 0: break
        n += 1
    return s
def kop(a, b):
    r = [0]*(N+1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i+j > N: break
                if y: r[i+j] += x*y
    return r
def teskari(a):  # a0=1 bo'lgan seriya teskarisi
    assert a[0] == 1
    r = [0]*(N+1); r[0] = 1
    for n in range(1, N+1):
        r[n] = -sum(a[k]*r[n-k] for k in range(1, n+1))
    return r
def darja(s, k):  # s^k (k manfiy bo'lishi mumkin)
    if k == 0:
        e = [0]*(N+1); e[0] = 1; return e
    b = s if k > 0 else teskari(s)
    r = b
    for _ in range(abs(k)-1): r = kop(r, b)
    return r

f14 = theta(1,4); f23 = theta(2,3); f312 = theta(3,12); f69 = theta(6,9); E5 = theta(5,10)
BAZA = [("f14",f14),("f23",f23),("f312",f312),("f69",f69),("E5",E5)]

# Monomiallar: darajalar -1..3, yig'indi 3
import itertools
monlar = []
for eks in itertools.product(range(-1,4), repeat=5):
    if sum(eks) != 3: continue
    if sum(1 for e in eks if e) > 3: continue  # eng ko'p 3 funksiya qatnashsin (kitob uslubi)
    s = [0]*(N+1); s[0] = 1
    for (nom,f),e in zip(BAZA,eks):
        if e: s = kop(s, darja(f,e))
    nomi = "*".join(f"{nom}^{e}" if e!=1 else nom for (nom,_),e in zip(BAZA,eks) if e)
    monlar.append((nomi, s))
print("monomiallar:", len(monlar))

# Hash: ±q^b * M (b=0..3) -> birinchi 80 koef tuple
H = {}
for nomi, s in monlar:
    for b in range(0, 4):
        for c in (1,-1):
            key = tuple(c*s[i-b] if i>=b else 0 for i in range(80))
            H.setdefault(key, f"{'+' if c>0 else '-'}q^{b}*{nomi}")

def qidir(T, nomT):
    print(f"\n### {nomT}")
    top = []
    for nomi, s in monlar:
        for a in range(0, 4):
            for c in (1,-1):
                R = [T[i] - c*(s[i-a] if i>=a else 0) for i in range(N+1)]
                key = tuple(R[:80])
                if key in H:
                    yech = f"{'+' if c>0 else '-'}q^{a}*{nomi}  {H[key]}"
                    if yech not in top:
                        top.append(yech); print("  YECHIM:", yech)
    if not top: print("  2-hadli topilmadi (bu bazada)")
    return top

# NAZORAT: (A1-B1)^3 = 8.2.9 = f23*f14^3/f312 + q^2*f312^3  (kitob)
A1 = theta(7,8); B1t = theta(2,13)
B1 = [0]*(N+1)
for i,x in enumerate(B1t):
    if i+1<=N: B1[i+1]=x
AB1m = [x-y for x,y in zip(A1,B1)]
T_naz = kop(kop(AB1m,AB1m),AB1m)
qidir(T_naz, "NAZORAT (A1-B1)^3 — kutiladi: +q^0*f23*f14^3*f312^-1 va +q^2*f312^3")

# NISHON 1: A1^3 - B1^3 (Entry 8.2.1 yo'qolgan izoh)
T1 = [x-y for x,y in zip(kop(kop(A1,A1),A1), kop(kop(B1,B1),B1))]
qidir(T1, "NISHON-1: A1^3-B1^3 (7,8)")

# NISHON 2: A2^3 + B2^3 (Entry 8.2.2 yo'qolgan izoh)
A2 = theta(4,11); B2t = theta(1,14)
B2 = [0]*(N+1)
for i,x in enumerate(B2t):
    if i+1<=N: B2[i+1]=x
T2 = [x+y for x,y in zip(kop(kop(A2,A2),A2), kop(kop(B2,B2),B2))]
qidir(T2, "NISHON-2: A2^3+B2^3 (4,11)")
