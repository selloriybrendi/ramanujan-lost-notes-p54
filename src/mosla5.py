# Buzuvchi-sinov: butun-q bazada 2-hadli, KENG koeffitsiyentlar (c=+-1..6), trivial istisno
import itertools
N = 200; K = 80
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
        e=[0]*(N+1); e[0]=1; return e
    b = s if k>0 else teskari(s)
    r = b
    for _ in range(abs(k)-1): r = kop(r, b)
    return r
f14=theta(1,4); f23=theta(2,3); f312=theta(3,12); f69=theta(6,9); E5=theta(5,10)
f114=theta(1,14); f213=theta(2,13); f411=theta(4,11); f78=theta(7,8)
E1=theta(1,2); E3=theta(3,6); E15=theta(15,30)
BAZA=[("f14",f14),("f23",f23),("f312",f312),("f69",f69),("E5",E5),
      ("f114",f114),("f213",f213),("f411",f411),("f78",f78),("E1",E1),("E3",E3),("E15",E15)]
monlar=[]
for r in (1,2,3):
    for tanl in itertools.combinations(range(len(BAZA)), r):
        for eks in itertools.product(range(-1,4), repeat=r):
            if sum(eks)!=3 or 0 in eks: continue
            s=[0]*(N+1); s[0]=1
            for idx,e in zip(tanl,eks): s=kop(s,darja(BAZA[idx][1],e))
            nomi="*".join(f"{BAZA[i][0]}^{e}" if e!=1 else BAZA[i][0] for i,e in zip(tanl,eks))
            monlar.append((nomi,s[:K]))
H={}
for nomi,s in monlar:
    for b in range(0,4):
        for c in (1,-1,2,-2,3,-3):
            key=tuple(c*s[i-b] if i>=b else 0 for i in range(K))
            H.setdefault(key,(c,b,nomi))
print("monomial:",len(monlar),"hash:",len(H))
A1=f78; B1=[0]*(N+1)
for i,x in enumerate(f213):
    if i+1<=N: B1[i+1]=x
A2=f411; B2=[0]*(N+1)
for i,x in enumerate(f114):
    if i+1<=N: B2[i+1]=x
T1=[x-y for x,y in zip(kop(kop(A1,A1),A1),kop(kop(B1,B1),B1))][:K]
T2=[x+y for x,y in zip(kop(kop(A2,A2),A2),kop(kop(B2,B2),B2))][:K]
def qidir(T,nomT,triv):
    top=[]
    for nomi,s in monlar:
        for a in range(0,4):
            for c in (1,-1,2,-2,3,-3,6,-6):
                key=tuple(T[i]-c*(s[i-a] if i>=a else 0) for i in range(K))
                if key in H:
                    c2,b2,n2=H[key]
                    y=tuple(sorted([f"{c:+d}*q^{a}*{nomi}",f"{c2:+d}*q^{b2}*{n2}"]))
                    if any(t in " ".join(y) for t in triv): continue
                    if y not in top: top.append(y); print(f"  {nomT}:"," ".join(y))
    if not top: print(f"  {nomT}: keng-c bilan ham 2-hadli chiqmadi")
qidir(T1,"N1",["f78^3","f213^3"])
qidir(T2,"N2",["f411^3","f114^3"])
