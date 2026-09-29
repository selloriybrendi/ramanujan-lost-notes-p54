# IJOD-4: yo'qolgan izohlar = Entry 8.2.3 (circular summation, O'SHA p.54) ning xususiy holi?
#  8.2.1: a=-q^3, b=-q^2 (ab=q^5):
#    A^3 - B^3 = q^2 f^3(-q^3,-q^12) + f(-q^2,-q^3) * { psi^3(q^5)/psi(q^15) + 3 q^5 psi^3(q^15)/psi(q^5) }
#  8.2.2: a=-q^6, b=-q^(-1) (ab=q^5):
#    A^3 + B^3 = q^(-1) [ f^3(-q^6,-q^9) - f(-q,-q^4) * { psi^3(q^5)/psi(q^15) + 3 q^5 psi^3(q^15)/psi(q^5) } ]
# Tekshiruv: bo'lishsiz, psi(q^5)psi(q^15) ga ko'paytirib, butun sonli qatorlarda.
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
K = N + 20
def theta(a, b):  # f(-q^a,-q^b)
    s=[0]*K; n=0
    while True:
        did=False
        for m in ((n,-n) if n else (0,)):
            e=a*m*(m+1)//2+b*m*(m-1)//2
            if 0<=e<K: s[e]+=(-1)**(m%2); did=True
        if not did and n>0: break
        n+=1
    return s
def psi(k):  # psi(q^k) = sum_{n>=0} q^{k n(n+1)/2}
    s=[0]*K; n=0
    while k*n*(n+1)//2 < K: s[k*n*(n+1)//2]+=1; n+=1
    return s
def kop(a,b):
    r=[0]*K
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:K-i]):
                if y: r[i+j]+=x*y
    return r
def kub(a): return kop(kop(a,a),a)
def sil(s,k): return [0]*k+s[:K-k]
def q(a,b): return [x+y for x,y in zip(a,b)]
def ay(a,b): return [x-y for x,y in zip(a,b)]
def sk(c,a): return [c*x for x in a]
P5, P15 = psi(5), psi(15)
# {psi^3(q^5)/psi(q^15) + 3q^5 psi^3(q^15)/psi(q^5)} * psi(q^5)psi(q^15) = psi^4(q^5) + 3 q^5 psi^4(q^15)
Kqavs = q(kop(kub(P5),P5), sil(sk(3,kop(kub(P15),P15)),5))
den = kop(P5,P15)
A1,B1 = theta(7,8), sil(theta(2,13),1)
A2,B2 = theta(4,11), sil(theta(1,14),1)
# 8.2.1: (A^3-B^3 - q^2 f312^3) * psi5 psi15 == f23 * Kqavs
L1 = kop(ay(ay(kub(A1),kub(B1)), sil(kub(theta(3,12)),2)), den)
R1 = kop(theta(2,3), Kqavs)
# 8.2.2: (q(A^3+B^3) - f69^3) * psi5 psi15 == - f14 * Kqavs
L2 = kop(ay(sil(q(kub(A2),kub(B2)),1), kub(theta(6,9))), den)
R2 = sk(-1, kop(theta(1,4), Kqavs))
for nom,L,R in (("8.2.1 lost: A^3-B^3 (circular summation)",L1,R1),("8.2.2 lost: A^3+B^3 (circular summation)",L2,R2)):
    ok = L[:N]==R[:N]
    d = None if ok else next(i for i in range(N) if L[i]!=R[i])
    print(f"{nom}: {'MOS ✓' if ok else 'FARQ ✗ q^%d: %d vs %d'%(d,L[d],R[d])}  ({N} had)")
