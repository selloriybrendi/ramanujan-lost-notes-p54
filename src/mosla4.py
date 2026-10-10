# q^(1/3)-panjara: hammasi q->q^3 dunyosida. Baza = 15-juftliklar (asl, "1/3" vakillari) + ularning 3x versiyalari.
# Nazorat: A1-B1 (q^3) = f(-q^2,-q^3) + q^2*f(-q^9,-q^36)  [8.2.7]
import itertools
K = 121  # koef oynasi
def theta(c4, c11):
    s = [0]*K; n = 0
    while True:
        did = False
        for m in (n, -n) if n else (0,):
            e = c4*m*(m+1)//2 + c11*m*(m-1)//2
            if 0 <= e < K: s[e] += (-1)**(m % 2); did = True
        if not did and n > 0: break
        n += 1
    return s
def kop(a, b):
    r = [0]*K
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i+j >= K: break
                if y: r[i+j] += x*y
    return r
def teskari(a):
    r = [0]*K; r[0] = 1
    for n in range(1, K): r[n] = -sum(a[k]*r[n-k] for k in range(1, n+1))
    return r
def darja(s, k):
    if k == 0:
        e=[0]*K; e[0]=1; return e
    b = s if k>0 else teskari(s)
    r = b
    for _ in range(abs(k)-1): r = kop(r, b)
    return r

juft = [(1,14),(2,13),(3,12),(4,11),(5,10),(6,9),(7,8),(1,4),(2,3),(1,2),(3,6),(15,30)]
hamma = sorted(set(juft + [(3*a,3*b) for a,b in juft]))
BAZA = [(f"F{a}_{b}", theta(a,b)) for a,b in hamma]
print("funksiyalar:", len(BAZA))

monlar = []
for r in (1,2,3):
    for tanl in itertools.combinations(range(len(BAZA)), r):
        for eks in itertools.product(range(-1,4), repeat=r):
            if sum(eks) != 3 or 0 in eks: continue
            s = [0]*K; s[0] = 1
            for idx,e in zip(tanl,eks): s = kop(s, darja(BAZA[idx][1], e))
            nomi = "*".join(f"{BAZA[i][0]}^{e}" if e!=1 else BAZA[i][0] for i,e in zip(tanl,eks))
            monlar.append((nomi, s))
print("monomiallar:", len(monlar))

import array
def kalit(v): return array.array('q', v).tobytes()
print("koeffitsientlar: c in {±1,±2,±3}; siljish 0..9; oyna K =", K)
H = {}
for nomi, s in monlar:
    for b in range(0, 10):
        for c in (1,-1,2,-2,3,-3):
            key = kalit([c*s[i-b] if i>=b else 0 for i in range(K)])
            H.setdefault(key, (c,b,nomi))
print("hash:", len(H))

def q3(s):  # s(q)->s(q^3), K oynada
    r = [0]*K
    for i,x in enumerate(s):
        if 3*i < K: r[3*i] = x
    return r

A1_3 = theta(21,24); B1_3 = [0]*K
t = theta(6,39)
for i,x in enumerate(t):
    if i+3 < K: B1_3[i+3] = x
A2_3 = theta(12,33); B2_3 = [0]*K
t = theta(3,42)
for i,x in enumerate(t):
    if i+3 < K: B2_3[i+3] = x

def qidir(T, nomT, trivialso):
    top = []
    for nomi, s in monlar:
        for a in range(0, 10):
            for c in (1,-1,2,-2,3,-3):
                key = kalit([T[i]-c*(s[i-a] if i>=a else 0) for i in range(K)])
                if key in H:
                    c2,b2,n2 = H[key]
                    y = tuple(sorted([f"{c:+d}*q^{a}*{nomi}", f"{c2:+d}*q^{b2}*{n2}"]))
                    if any(t in " ".join(y) for t in trivialso): continue
                    if y not in top:
                        top.append(y); print(f"  {nomT} YECHIM:", " ".join(y))
    if not top: print(f"  {nomT}: 2-hadli chiqmadi")

# 1-darajali NAZORAT: A1-B1 (bu yerda monomial darajasi 1 kerak — alohida kichik qidiruv)
mon1 = []
for idx,(nomi,s) in enumerate(BAZA):
    mon1.append((nomi, s))
H1 = {}
for nomi,s in mon1:
    for b in range(0,10):
        for c in (1,-1):
            key = tuple(c*s[i-b] if i>=b else 0 for i in range(K))
            H1.setdefault(key,(c,b,nomi))
T_naz = [x-y for x,y in zip(A1_3,B1_3)]
naz_ok = False
for nomi,s in mon1:
    for a in range(0,10):
        for c in (1,-1):
            key = tuple(T_naz[i]-c*(s[i-a] if i>=a else 0) for i in range(K))
            if key in H1:
                c2,b2,n2 = H1[key]
                print("NAZORAT (A1-B1)@q^3:", f"{c:+d}*q^{a}*{nomi}  {c2:+d}*q^{b2}*{n2}")
                naz_ok = True
print("nazorat topildi:", naz_ok)

T1 = [x-y for x,y in zip(kop(kop(A1_3,A1_3),A1_3), kop(kop(B1_3,B1_3),B1_3))]
T2 = [x+y for x,y in zip(kop(kop(A2_3,A2_3),A2_3), kop(kop(B2_3,B2_3),B2_3))]
qidir(T1, "N1 (7,8) A^3-B^3 @q^(1/3)", ["F21_24^3","F6_39^3"])
qidir(T2, "N2 (4,11) A^3+B^3 @q^(1/3)", ["F12_33^3","F3_42^3"])
