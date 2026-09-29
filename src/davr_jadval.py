# davr_jadval.py — maqola §3.1 jadvali uchun o'lchov:
# 4 kvadratik kofaktorning Euler-ko'rsatkichlari (n<=380), davr 1..60 sinovi.
# Chiqish (2026-09-29 yurgizilgan):
#   8.2.1 A^2-AB+B^2 (yozilgan 8.2.8 omili): davr=15
#   8.2.1 A^2+AB+B^2 (yo'qolgan A^3-B^3 omili): davr=YO'Q (1..60, n<=380)
#   8.2.2 A^2+AB+B^2 (yozilgan 8.2.12 omili): davr=15
#   8.2.2 A^2-AB+B^2 (yo'qolgan A^3+B^3 omili): davr=YO'Q (1..60, n<=380)
N = 400
def theta(a,b):
    s=[0]*(N+1)
    n=0
    while True:
        did=False
        for m in ((n,-n) if n else (0,)):
            e=a*m*(m+1)//2+b*m*(m-1)//2
            if 0<=e<=N:
                s[e]+=(-1)**(m%2)
                did=True
        if not did and n>0: break
        n+=1
    return s
def kop(a,b):
    r=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]):
                if y: r[i+j]+=x*y
    return r
def euler(s):
    v=next(i for i,x in enumerate(s) if x)
    u=s[v]
    assert abs(u)==1
    c=[x*u for x in s[v:]] if u==-1 else s[v:]
    M=len(c)-1
    num=[i*c[i] for i in range(len(c))]
    L=[0]*(M+1)
    for m in range(M+1):
        t=num[m]
        for k in range(1,m+1):
            t-=c[k]*L[m-k] if k<len(c) else 0
        L[m]=t//c[0]
    e=[0]*(M+1)
    for m in range(1,M+1):
        t=-L[m]
        for d in range(1,m):
            if m%d==0: t-=d*e[d]
        assert t%m==0
        e[m]=t//m
    return e[1:]
def davr_top(e, nmax=380, dmax=60):
    e=e[:nmax]
    for d in range(1,dmax+1):
        if all(e[i]==e[i+d] for i in range(len(e)-d-20)):
            return d
    return None
def sur(t):
    B=[0]*(N+1)
    for i,x in enumerate(t):
        if i+1<=N: B[i+1]=x
    return B
A1=theta(7,8); B1=sur(theta(2,13))
A2=theta(4,11); B2=sur(theta(1,14))
def kvad(A,B,sgn):
    return [a+sgn*m+b for a,m,b in zip(kop(A,A),kop(A,B),kop(B,B))]
if __name__=="__main__":
    for nom,A,B,sgn in [("8.2.1 A^2-AB+B^2 (yozilgan 8.2.8 omili)",A1,B1,-1),
                        ("8.2.1 A^2+AB+B^2 (yo'qolgan A^3-B^3 omili)",A1,B1,+1),
                        ("8.2.2 A^2+AB+B^2 (yozilgan 8.2.12 omili)",A2,B2,+1),
                        ("8.2.2 A^2-AB+B^2 (yo'qolgan A^3+B^3 omili)",A2,B2,-1)]:
        d=davr_top(euler(kvad(A,B,sgn)))
        print(f"{nom}: davr={'YO`Q (1..60 sinaldi, n<=380)' if d is None else d}")
