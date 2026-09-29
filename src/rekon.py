# Rekonstruksiya-tasdiqlash (q^3 dunyoda, 120 had):
# N1: A^3-B^3 = [f23*f14^3/f312 + q^2*f312^3]@q^3  +  3*AB*[F2_3 + q^2*F9_36]
# N2: A^3+B^3 = [ -1/q {f14*f23^3/f69 - f69^3} ]@q^3  -  3*AB*[ -1/q^(1/3){F1_4... } ]  (8.2.13 - 3AB*8.2.11)
K = 121
def theta(a, b):
    s = [0]*K; n = 0
    while True:
        did = False
        for m in (n, -n) if n else (0,):
            e = a*m*(m+1)//2 + b*m*(m-1)//2
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
def sil(s, k):  # q^k ga ko'paytirish
    r = [0]*K
    for i,x in enumerate(s):
        if i+k < K and i+k >= 0: r[i+k] = x
    return r
def ayr(a,b): return [x-y for x,y in zip(a,b)]
def qosh(a,b): return [x+y for x,y in zip(a,b)]

# q^3-dunyoda funksiyalar
F14=theta(3,12); F23=theta(6,9); F312=theta(9,36); F69=theta(18,27)
F2_3=theta(2,3); F9_36=theta(9,36)  # 8.2.7 hadlari (q^(1/3) vakillari)
A1=theta(21,24); B1=sil(theta(6,39),3)
A2=theta(12,33); B2=sil(theta(3,42),3)
AB1=kop(A1,B1); AB2=kop(A2,B2)

# N1
T1 = ayr(kop(kop(A1,A1),A1), kop(kop(B1,B1),B1))
had829 = qosh(kop(kop(F23,kop(kop(F14,F14),F14)),teskari(F312)), sil(kop(kop(F312,F312),F312),6))  # q^2@q^3=q^6
AmB = qosh(F2_3, sil(F9_36,2))
R1 = qosh(had829, kop([3*x for x in AB1], AmB))
print("N1: A^3-B^3 == (8.2.9) + 3AB*(8.2.7):", "MOS ✓" if T1==R1 else "FARQ ✗")

# N2: 8.2.13: (A+B)^3 = -1/q { f14 f23^3 / f69 - f69^3 }  (butun-q tilda) -> q^3: -q^{-3}{F14 F23^3/F69 - F69^3}
# 8.2.11: A+B = -q^{-1/3} { f(-q^{1/3},-q^{4/3}) - f69 } -> q^3: -q^{-1}{ theta(1,4) - F69 }
ApB3 = [ -x for x in sil(ayr(kop(kop(F14,kop(kop(F23,F23),F23)),teskari(F69)), kop(kop(F69,F69),F69)), -3) ]
ApB  = [ -x for x in sil(ayr(theta(1,4), F69), -1) ]
T2 = qosh(kop(kop(A2,A2),A2), kop(kop(B2,B2),B2))
R2 = ayr(ApB3, kop([3*x for x in AB2], ApB))
print("N2: A^3+B^3 == (8.2.13) - 3AB*(8.2.11):", "MOS ✓" if T2==R2 else "FARQ ✗")
if T2!=R2:
    d=next(i for i in range(K) if T2[i]!=R2[i]); print("farq q^",d,T2[d],R2[d])
# Sanity: (A+B) seriya to'g'riligini alohida
print("8.2.11 sanity (A2+B2==ApB):", "MOS ✓" if qosh(A2,B2)==ApB else "FARQ ✗")
print("8.2.13 sanity ((A2+B2)^3==ApB3):", "MOS ✓" if kop(kop(qosh(A2,B2),qosh(A2,B2)),qosh(A2,B2))==ApB3 else "FARQ ✗")
print("8.2.7 sanity (A1-B1==AmB):", "MOS ✓" if ayr(A1,B1)==AmB else "FARQ ✗")
