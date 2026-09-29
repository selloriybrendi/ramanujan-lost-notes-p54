# sonli_tekshir.py — maqoladagi formulalarni YOZILGAN ko'rinishida (bo'lish, q^(1/3)) 50 xonali aniqlikda hisoblaydi.
# verifier.py (butun-son qator, ko'paytirib solishtirish) dan MUSTAQIL usul.
from mpmath import mp, mpf, nsum, inf
mp.dps = 50
def f(a, b):  # Ramanujan f(a,b) = sum a^{n(n+1)/2} b^{n(n-1)/2}
    s = mpf(0); n = 0
    while True:
        t1 = a**(n*(n+1)//2) * b**(n*(n-1)//2)
        t2 = a**((-n)*(-n+1)//2) * b**((-n)*(-n-1)//2) if n else 0
        s += t1 + t2
        if abs(t1) + abs(t2) < mpf(10)**(-60) and n > 3: break
        n += 1
    return s
def poch(a, q):  # (a;q)_inf
    p = mpf(1); k = 0
    while True:
        t = a * q**k
        if abs(t) < mpf(10)**(-60): break
        p *= (1 - t); k += 1
    return p
E = lambda x: f(-x, -x*x)  # f(-x)
def a_(x):  # Borwein a(x) = sum_{m,n} x^{m^2+mn+n^2}
    s = mpf(0); M = 1
    while x**(M*M*3//4) > mpf(10)**-60: M += 1
    for m in range(-M, M+1):
        for n in range(-M, M+1): s += x**(m*m+m*n+n*n)
    return s
worst = mpf(0)
for q in [mpf('0.05'), mpf('0.2'), mpf('0.37'), mpf('0.55'), mpf('0.8')]:
    r = q**(mpf(1)/3)
    A1, B1 = f(-q**7, -q**8), q*f(-q**2, -q**13)
    A2, B2 = f(-q**4, -q**11), q*f(-q, -q**14)
    f14, f23, f312, f69 = f(-q, -q**4), f(-q**2, -q**3), f(-q**3, -q**12), f(-q**6, -q**9)
    t = {
     '8.2.6':  (A1+B1, f23/f14*E(q**5)),
     '8.2.7':  (A1-B1, f(-r**2, -q) + r**2*f312),
     '8.2.8':  (A1**3+B1**3, f69/f312*E(q**5)**3),
     '8.2.9':  ((A1-B1)**3, f23*f14**3/f312 + q**2*f312**3),
     'T1':     (A1**3-B1**3, f23*f14**3/f312 + q**2*f312**3 + 3*A1*B1*(f(-r**2, -q) + r**2*f312)),
     'T1-AB':  (A1*B1, q*poch(q**2,q**15)*poch(q**7,q**15)*poch(q**8,q**15)*poch(q**13,q**15)*poch(q**15,q**15)**2),
     '8.2.10': (A2-B2, f14/f23*E(q**5)),
     '8.2.11': (A2+B2, -(f(-r, -r**4) - f69)/r),
     '8.2.12': (A2**3-B2**3, f312/f69*E(q**5)**3),
     '8.2.13': ((A2+B2)**3, -(f14*f23**3/f69 - f69**3)/q),
     'T2':     (A2**3+B2**3, -(f14*f23**3/f69 - f69**3)/q + 3*A2*B2/r*(f(-r, -r**4) - f69)),
     'T3':     (A1**3-B1**3, q**2*f312**3 + f23*a_(q**5)),
     'T4':     (A2**3+B2**3, (f69**3 - f14*a_(q**5))/q),
     'T2-AB':  (A2*B2, q*poch(q,q**15)*poch(q**4,q**15)*poch(q**11,q**15)*poch(q**14,q**15)*poch(q**15,q**15)**2),
    }
    for k, (L, R) in t.items():
        e = abs(L-R)/max(abs(L), mpf(1)); worst = max(worst, e)
    print(f"q={q}: 14 formula, max nisbiy farq =", mp.nstr(max(abs(L-R)/max(abs(L),1) for L,R in t.values()), 3))
print("JAMI eng katta farq:", mp.nstr(worst, 3), "->", "MOS" if worst < mpf(10)**(-40) else "FARQ")
