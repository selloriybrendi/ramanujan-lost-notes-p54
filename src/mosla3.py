# Nazariy sinov: T ± 3*q^a*(aralash uch-mahsulot) ni 2-hadli hashda qidirish
import itertools
N = 260; K = 80
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
T1=[x-y for x,y in zip(kop(kop(A1,A1),A1),kop(kop(B1,B1),B1))]
T2=[x+y for x,y in zip(kop(kop(A2,A2),A2),kop(kop(B2,B2),B2))]

# aralash hadlar: AB*f312, AB*f69 (har juftlik uchun), q-siljishlar bilan
AB1=kop(A1,B1); AB2=kop(A2,B2)
aralash1=[("AB1*f312",kop(AB1,f312)[:K]),("AB1*f69",kop(AB1,f69)[:K]),("AB1*E5",kop(AB1,E5)[:K])]
aralash2=[("AB2*f312",kop(AB2,f312)[:K]),("AB2*f69",kop(AB2,f69)[:K]),("AB2*E5",kop(AB2,E5)[:K])]

def uchinchi(T,arlar,nomT,trivial):
    Tk=T[:K]; top=[]
    for nomi,s in arlar:
        for a in range(0,3):
            for c in (3,-3):
                R=[Tk[i]-c*(s[i-a] if i>=a else 0) for i in range(K)]
                # R ni 2-hadli qidir
                for n2,s2 in monlar:
                    for a2 in range(0,4):
                        for c2 in (1,-1,2,-2,3,-3):
                            key=tuple(R[i]-c2*(s2[i-a2] if i>=a2 else 0) for i in range(K))
                            if key in H:
                                c3,b3,n3=H[key]
                                yech=(f"{c:+d}*q^{a}*{nomi}",)+tuple(sorted([f"{c2:+d}*q^{a2}*{n2}",f"{c3:+d}*q^{b3}*{n3}"]))
                                if trivial(yech): continue
                                if yech not in top:
                                    top.append(yech); print(f"  {nomT} YECHIM:", " ".join(yech))
    if not top: print(f"  {nomT}: topilmadi")
uch1=lambda y: ("f78^3" in " ".join(y) or "f213^3" in " ".join(y))
uch2=lambda y: ("f411^3" in " ".join(y) or "f114^3" in " ".join(y))
print("== NISHON-1 (7,8) A^3-B^3, 3-hadli (nazariy AB-had bilan):")
uchinchi(T1,aralash1,"N1",uch1)
print("== NISHON-2 (4,11) A^3+B^3:")
uchinchi(T2,aralash2,"N2",uch2)
