# IJOD-5: yo'qolgan izohlar Borwein a(q) orqali (bo'lishsiz, to'g'ridan a(q)=sum q^{m^2+mn+n^2}):
#   (N1)  A^3 - B^3 = q^2 f^3(-q^3,-q^12) + f(-q^2,-q^3) a(q^5)          [A=f(-q^7,-q^8), B=q f(-q^2,-q^13)]
#   (N2)  q(A^3 + B^3) = f^3(-q^6,-q^9) - f(-q,-q^4) a(q^5)              [A=f(-q^4,-q^11), B=q f(-q,-q^14)]
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
K = N + 5
def theta(a,b):
    s=[0]*K; n=0
    while True:
        did=False
        for m in ((n,-n) if n else (0,)):
            e=a*m*(m+1)//2+b*m*(m-1)//2
            if 0<=e<K: s[e]+=(-1)**(m%2); did=True
        if not did and n>0: break
        n+=1
    return s
def a_borwein(k):  # a(q^k)
    s=[0]*K; M=int((4*K/(3*k))**0.5)+3  # m^2+mn+n^2 >= (3/4)max^2
    for m in range(-M,M+1):
        for n in range(-M,M+1):
            e=k*(m*m+m*n+n*n)
            if e<K: s[e]+=1
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
a5=a_borwein(5)
A1,B1=theta(7,8),sil(theta(2,13),1); A2,B2=theta(4,11),sil(theta(1,14),1)
L1=[x-y for x,y in zip(kub(A1),kub(B1))]; R1=[x+y for x,y in zip(sil(kub(theta(3,12)),2), kop(theta(2,3),a5))]
L2=sil([x+y for x,y in zip(kub(A2),kub(B2))],1); R2=[x-y for x,y in zip(kub(theta(6,9)), kop(theta(1,4),a5))]
for nom,L,R in (("N1: A^3-B^3 = q^2 f^3(-q^3,-q^12) + f(-q^2,-q^3) a(q^5)",L1,R1),("N2: q(A^3+B^3) = f^3(-q^6,-q^9) - f(-q,-q^4) a(q^5)",L2,R2)):
    ok=L[:N]==R[:N]; d=None if ok else next(i for i in range(N) if L[i]!=R[i])
    print(f"{nom}: {'MOS ✓' if ok else 'FARQ q^%d'%d}  ({N} had)")
# manfiy nazorat: a(q^5) o'rniga a(q^5) ning bitta koeffitsiyenti buzilsa
a5b=a5[:]; a5b[10]+=1
R1b=[x+y for x,y in zip(sil(kub(theta(3,12)),2), kop(theta(2,3),a5b))]
print("nazorat (buzilgan a):", "FARQ ✓ (test sezgir)" if L1[:N]!=R1b[:N] else "MOS ✗ (test ko'r!)")
