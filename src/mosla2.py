# Kengaytirilgan baza + 2-hadli, keyin 3-hadli qidiruv (80-koef hash fazosida, so'ng to'liq tasdiqlash)
import itertools
N = 260; K = 80
MON = {}  # nomi -> to'liq N-koeffitsientli qator (tasdiq uchun)
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
def teskari(a):
    r = [0]*(N+1); r[0] = 1
    for n in range(1, N+1): r[n] = -sum(a[k]*r[n-k] for k in range(1, n+1))
    return r
def darja(s, k):
    if k == 0:
        e = [0]*(N+1); e[0] = 1; return e
    b = s if k > 0 else teskari(s)
    r = b
    for _ in range(abs(k)-1): r = kop(r, b)
    return r

BAZA = [("f14",theta(1,4)),("f23",theta(2,3)),("f312",theta(3,12)),("f69",theta(6,9)),
        ("E5",theta(5,10)),("f114",theta(1,14)),("f213",theta(2,13)),("f411",theta(4,11)),
        ("f78",theta(7,8)),("f510",theta(5,10+5)),("E1",theta(1,2)),("E3",theta(3,6)),("E15",theta(15,30))]
# f510 xato bo'lmasin: f(-q^5,-q^10)=theta(5,10) — lekin E5 ham theta(5,10)! To'g'irlash:
# f(-q^5) Euler = f(-q^5,-q^10) AYNAN bir xil narsa (f(-x)=f(-x,-x^2)). Demak E5 ni bitta qoldiramiz.
BAZA = [("f14",theta(1,4)),("f23",theta(2,3)),("f312",theta(3,12)),("f69",theta(6,9)),
        ("E5",theta(5,10)),("f114",theta(1,14)),("f213",theta(2,13)),("f411",theta(4,11)),
        ("f78",theta(7,8)),("E1",theta(1,2)),("E3",theta(3,6)),("E15",theta(15,30))]

monlar = []
for tanl in itertools.combinations(range(len(BAZA)), 3):
    for eks in itertools.product(range(-1,4), repeat=3):
        if sum(eks) != 3 or 0 in eks: continue
        s = [0]*(N+1); s[0] = 1
        for idx, e in zip(tanl, eks): s = kop(s, darja(BAZA[idx][1], e))
        nomi = "*".join(f"{BAZA[i][0]}^{e}" if e!=1 else BAZA[i][0] for i,e in zip(tanl,eks))
        monlar.append((nomi, s[:K])); MON[nomi] = s
# 1 va 2 funksiyali monomiallar ham
for tanl in list(itertools.combinations(range(len(BAZA)),1))+list(itertools.combinations(range(len(BAZA)),2)):
    for eks in itertools.product(range(-1,4), repeat=len(tanl)):
        if sum(eks) != 3 or 0 in eks: continue
        s = [0]*(N+1); s[0] = 1
        for idx, e in zip(tanl, eks): s = kop(s, darja(BAZA[idx][1], e))
        nomi = "*".join(f"{BAZA[i][0]}^{e}" if e!=1 else BAZA[i][0] for i,e in zip(tanl,eks))
        monlar.append((nomi, s[:K])); MON[nomi] = s
print("monomiallar:", len(monlar))

print("koeffitsientlar: c in {±1,±2,±3}; siljish 0..3; hash K =", K, "; tasdiq N =", N)
H = {}
for nomi, s in monlar:
    for b in range(0, 4):
        for c in (1,-1,2,-2,3,-3):
            key = tuple(c*s[i-b] if i>=b else 0 for i in range(K))
            H.setdefault(key, (c,b,nomi))

def qidir2(T, nomT):
    Tk = T[:K]; top = []
    for nomi, s in monlar:
        for a in range(0, 4):
            for c in (1,-1,2,-2,3,-3):
                key = tuple(Tk[i] - c*(s[i-a] if i>=a else 0) for i in range(K))
                if key in H:
                    c2,b2,n2 = H[key]
                    # to'liq tasdiq: barcha N koeffitsientda (hash faqat K=80 ni ko'radi)
                    s1 = MON[nomi]; s2 = MON[n2]; tas = all(T[i] == c*(s1[i-a] if i>=a else 0) + c2*(s2[i-b2] if i>=b2 else 0) for i in range(N+1))
                    yech = tuple(sorted([f"{c:+d}*q^{a}*{nomi}", f"{c2:+d}*q^{b2}*{n2}"])) + (("[260-koeff tasdiq: MOS]" if tas else "[260-koeff tasdiq: FARQ — rad]"),)
                    if yech not in top: top.append(yech)
    print(f"### {nomT}: {len(top)} ta 2-hadli nomzod")
    for y in top[:6]: print("   ", " ".join(y))
    return top

A1 = theta(7,8); B1t = theta(2,13); B1=[0]*(N+1)
for i,x in enumerate(B1t):
    if i+1<=N: B1[i+1]=x
A2 = theta(4,11); B2t = theta(1,14); B2=[0]*(N+1)
for i,x in enumerate(B2t):
    if i+1<=N: B2[i+1]=x
T1 = [x-y for x,y in zip(kop(kop(A1,A1),A1), kop(kop(B1,B1),B1))]
T2 = [x+y for x,y in zip(kop(kop(A2,A2),A2), kop(kop(B2,B2),B2))]
# NAZORAT (12-bazada, xuddi shu qidiruv): (A1-B1)^3 = (8.2.9) qayta topilishi kerak
T0 = [x-y for x,y in zip(A1,B1)]; T0 = kop(kop(T0,T0),T0)
t0 = qidir2(T0, "NAZORAT: (A1-B1)^3 (7,8) — kutiladi (8.2.9): f23*f14^3/f312 + q^2*f312^3")
t1 = qidir2(T1, "NISHON-1: A1^3-B1^3 (7,8)")
t2 = qidir2(T2, "NISHON-2: A2^3+B2^3 (4,11)")
