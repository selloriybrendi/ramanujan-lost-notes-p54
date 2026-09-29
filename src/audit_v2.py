# audit_v2.py — maqola v2 dagi HAR formula va isbotning HAR qadami, 50 xona (mpmath), bir nechta q va (a,b) da.
from mpmath import mp, mpf, cbrt
mp.dps = 50
EPS = mpf(10)**-40
def f(a, b):
    s = mpf(0); n = 0
    while True:
        t1 = a**(n*(n+1)//2) * b**(n*(n-1)//2)
        t2 = a**((-n)*(-n+1)//2) * b**((-n)*(-n-1)//2) if n else 0
        s += t1 + t2
        if abs(t1) + abs(t2) < mpf(10)**-60 and n > 3: break
        n += 1
    return s
E = lambda x: f(-x, -x*x)
def a_(x):
    s = mpf(0); M = 1
    while x**(M*M*3//4) > mpf(10)**-60: M += 1
    for m in range(-M, M+1):
        for n in range(-M, M+1): s += x**(m*m+m*n+n*n)
    return s
def poch(a, q):
    p = mpf(1); k = 0
    while abs(a*q**k) > mpf(10)**-60: p *= 1 - a*q**k; k += 1
    return p
nat = []
def sina(nom, L, R):
    e = abs(L-R)/max(abs(L), 1); ok = e < EPS; nat.append(ok)
    print(f"  {'MOS ✓' if ok else 'FARQ ✗'}  {mp.nstr(e,3):>9}  {nom}")
print("== (eq:circ) Entry 8.2.3 — ixtiyoriy (a,b), |ab|<1:")
for a, b in [(mpf('0.3'), mpf('-0.5')), (mpf('-0.2'), mpf('0.7')), (mpf('0.45'), mpf('0.6'))]:
    q = a*b
    L = f(a*b*b, a*a*b)**3 + a*f(b, a**3*b*b)**3 + b*f(a, a*a*b**3)**3
    R = f(a, b) * cbrt(E(q)**9/E(q**3)**3 + 27*q*E(q**3)**9/E(q)**3)
    sina(f"8.2.3 a={a}, b={b}", L, R)
print("== Borwein: a^3=b^3+c^3, brace = a(q):")
for q in [mpf('0.1'), mpf('0.4'), mpf('0.7')]:
    bq = E(q)**3/E(q**3); cq = 3*q**(mpf(1)/3)*E(q**3)**3/E(q)
    sina(f"a^3 = b^3+c^3  q={q}", a_(q)**3, bq**3 + cq**3)
    sina(f"brace = a(q)   q={q}", cbrt(E(q)**9/E(q**3)**3 + 27*q*E(q**3)**9/E(q)**3), a_(q))
for q in [mpf('0.2'), mpf('0.5'), mpf('0.8')]:
    A1, B1 = f(-q**7, -q**8), q*f(-q**2, -q**13)
    A2, B2 = f(-q**4, -q**11), q*f(-q, -q**14)
    f14, f23, f312, f69 = f(-q, -q**4), f(-q**2, -q**3), f(-q**3, -q**12), f(-q**6, -q**9)
    print(f"== q={q}: Theorem 1 isbot qadamlari (a=-q^3, b=-q^2)")
    a, b = -q**3, -q**2
    sina("ab = q^5", a*b, q**5)
    sina("f(ab^2,a^2b) = A", f(a*b*b, a*a*b), A1)
    sina("a f^3(b,a^3b^2) = -B^3", a*f(b, a**3*b*b)**3, -B1**3)
    sina("b f^3(a,a^2b^3) = -q^2 f^3(-q^3,-q^12)", b*f(a, a*a*b**3)**3, -q**2*f312**3)
    sina("f(a,b) = f(-q^2,-q^3)", f(a, b), f23)
    sina("THEOREM 1", A1**3-B1**3, q**2*f312**3 + f23*a_(q**5))
    print(f"== q={q}: Theorem 2 isbot qadamlari (a=-q^6, b=-q^-1)")
    a, b = -q**6, -1/q
    sina("ab = q^5", a*b, q**5)
    sina("f(ab^2,a^2b) = A", f(a*b*b, a*a*b), A2)
    sina("f(-q^-1,-q^16) = -q^-1 f(-q,-q^14)", f(-1/q, -q**16), -f(-q, -q**14)/q)
    sina("f(-q^-1,-q^6) = -q^-1 f(-q,-q^4)", f(-1/q, -q**6), -f14/q)
    sina("a f^3(b,a^3b^2) = B^3", a*f(b, a**3*b*b)**3, B2**3)
    sina("b f^3(a,a^2b^3) = -q^-1 f^3(-q^6,-q^9)", b*f(a, a*a*b**3)**3, -f69**3/q)
    sina("f(a,b) = -q^-1 f(-q,-q^4)", f(a, b), -f14/q)
    sina("THEOREM 2", A2**3+B2**3, (f69**3 - f14*a_(q**5))/q)
    print(f"== q={q}: by-product va JTP ko'paytmalar")
    sina("3AB(A-B) = f23{a(q^5) - f14^3/f312}", 3*A1*B1*(A1-B1), f23*(a_(q**5) - f14**3/f312))
    sina("3AB(A+B) = q^-1 f14{a(q^5) - f23^3/f69}  [8.2.2]", 3*A2*B2*(A2+B2), f14*(a_(q**5) - f23**3/f69)/q)
    sina("NAZORAT (buzilgan: Theorem 2 da +f14 a)", A2**3+B2**3, (f69**3 + f14*a_(q**5))/q)
    sina("AB (8.2.1) = q(q2,q7,q8,q13;q15)(q15;q15)^2", A1*B1, q*poch(q**2,q**15)*poch(q**7,q**15)*poch(q**8,q**15)*poch(q**13,q**15)*poch(q**15,q**15)**2)
    sina("AB (8.2.2) = q(q,q4,q11,q14;q15)(q15;q15)^2", A2*B2, q*poch(q,q**15)*poch(q**4,q**15)*poch(q**11,q**15)*poch(q**14,q**15)*poch(q**15,q**15)**2)
kutilgan_farq = 3  # har q da bitta NAZORAT FARQ bo'lishi SHART
farq = len(nat) - sum(nat)
print(f"\nAUDIT: {sum(nat)}/{len(nat)} MOS · FARQ {farq} (kutilgan {kutilgan_farq} = manfiy nazoratlar) · {'TO\'G\'RI' if farq == kutilgan_farq else 'XATO'}")
