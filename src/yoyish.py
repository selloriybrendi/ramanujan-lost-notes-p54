# A=f(-q^4,-q^11), B=q*f(-q,-q^14) — q-yoyish, N gacha; Euler-ko'rsatkichlar davriyligi
N = 400
def theta(c4, c11, tesk=False):
    # f(a,b)=sum a^{n(n+1)/2} b^{n(n-1)/2}; bu yerda a=-q^{c4}, b=-q^{c11} -> koef (-1)^n, daraja (c4+c11)n^2/2 + (c4-c11)n/2... aniq: c4*n(n+1)/2 + c11*n(n-1)/2
    s = [0]*(N+1)
    n = 0
    while True:
        did = False
        for m in (n, -n) if n else (0,):
            e = c4*m*(m+1)//2 + c11*m*(m-1)//2
            if 0 <= e <= N:
                s[e] += (-1)**(m % 2)
                did = True
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

def euler_korsatkich(s):
    # s = q^v * u * prod (1-q^n)^{e_n}; qaytaradi (v,u,[e_1..])
    v = next(i for i, x in enumerate(s) if x)
    u = s[v]
    c = [x//u for x in s[v:]]  # normalize, boshi 1
    assert c[0] == 1 and all(x*u == y for x, y in zip(c, s[v:]))
    e = []
    for n in range(1, len(c)):
        en = -c[n]  # (1-q^n)^{-e} ga bo'lish: log usuli o'rniga iterativ
        e.append(en)
        # c ni (1-q^n)^{en} ga bo'lish = (1-q^n)^{-en} ga ko'paytirish
        # (1-q^n)^{-en} = sum C(en+k-1,k) q^{nk}
        yangi = c[:]
        from math import comb
        for k in range(1, len(c)//n + 1):
            coef = comb(en+k-1, k) if en > 0 else ((-1)**k * comb(-en, k) if -en >= k else 0)
            if coef == 0 and en < 0: continue
            for i in range(0, len(c)-n*k):
                if c[i]: yangi[i+n*k] += c[i]*coef
        # to'g'ri usul: ketma-ket bo'lish
        pass
    return v, u, e

# ISHONCHLI usul: log-hosila bilan: agar F = prod (1-q^n)^{e_n} bo'lsa, -q F'/F = sum e_n * n*q^n/(1-q^n) = sum_m (sum_{d|m} d*e_d) q^m
def log_korsatkich(s):
    v = next(i for i, x in enumerate(s) if x)
    u = s[v]
    c = s[v:]
    # L_m = koeff of q^m in -q*d/dq log(F/ (u q^v)) ... hisob: F=c (c0=u). q F' / F: F' koeff (i)*c[i]
    M = len(c)-1
    # A_m: q F'/F = (sum i*c_i q^i)/ (sum c_i q^i)
    num = [i*c[i] for i in range(len(c))]
    L = [0]*(M+1)
    inv0 = 1  # c0 = u, ishlaymiz ratsional emas — c0 ga bo'lamiz (butun bo'lmasligi mumkin, lekin u=±1 kutiladi)
    assert abs(u) == 1
    cc = [x*u for x in c] if u == -1 else c[:]  # endi cc0=1... u=-1 bo'lsa ham nisbat o'zgarmaydi aslida
    # L = num/den (formal): L_m = num_m - sum_{k=1}^{m} cc_k L_{m-k}
    for m in range(M+1):
        t = num[m] if u == 1 else -num[m]*0 + num[m]* (1 if u==1 else 1)
    # soddalik: den=c, num=[i*c_i]; L_m = (num_m - sum_{k=1}^m c_k L_{m-k}/?) — bo'lish c0 ga
    L = [0]*(M+1)
    for m in range(M+1):
        t = num[m]
        for k in range(1, m+1):
            t -= c[k]*L[m-k] if k < len(c) else 0
        assert t % c[0] == 0
        L[m] = t // c[0]
    # endi L_m = -(sum_{d|m} d e_d)  chunki log(1-q^n) hosilasi manfiy
    S = [0]*(M+1)
    e = [0]*(M+1)
    for m in range(1, M+1):
        t = -L[m]
        for d in range(1, m):
            if m % d == 0: t -= d*e[d]
        assert t % m == 0, (m, t)
        e[m] = t//m
    return v, u, e[1:]

def f_seriya(nom, s, moslama=""):
    v, u, e = log_korsatkich(s)
    print(f"== {nom}: v={v} u={u}")
    print("e_1..e_45:", e[:45])
    # davriylik testi: davr 15 va 30
    for davr in (15, 30):
        ok = all(e[i] == e[i+davr] for i in range(len(e)-davr-5))
        print(f"   davr {davr}: {'DAVRIY ✓' if ok else 'yo`q'}")
    return e

A = theta(4, 11)
Bq = theta(1, 14)
B = [0]*(N+1)
for i, x in enumerate(Bq):
    if i+1 <= N: B[i+1] = x

# NAZORAT: A-B ma'lum mahsulotga ega (Entry 8.2.10) -> davriy chiqishi SHART
AB1 = [x-y for x, y in zip(A, B)]
AB2 = [x+y for x, y in zip(A, B)]
f_seriya("A-B (nazorat, ma'lum ayniyat)", AB1)
f_seriya("A+B (nazorat, ma'lum ayniyat)", AB2)

A3 = kop(kop(A, A), A)
B3 = kop(kop(B, B), B)
A3B3 = [x-y for x, y in zip(A3, B3)]
print("\nA^3-B^3 birinchi 30 koef:", A3B3[:30])
e = f_seriya("A^3-B^3 (NISHON)", A3B3)
