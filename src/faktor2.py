# (7,8): A^2+AB+B^2 davriymi?  (4,11): A^2-AB+B^2 davriymi?
N = 400
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
def log_e(s):
    v = next(i for i, x in enumerate(s) if x); c = s[v:]
    num = [i*c[i] for i in range(len(c))]; M = len(c)-1; L=[0]*(M+1)
    for m in range(M+1):
        t = num[m]
        for k in range(1, m+1): t -= c[k]*L[m-k] if k < len(c) else 0
        L[m] = t//c[0]
    e=[0]*(M+1)
    for m in range(1, M+1):
        t = -L[m]
        for d in range(1, m):
            if m % d == 0: t -= d*e[d]
        e[m] = t//m
    return v, s[v], e[1:]
def sinov(nom, s):
    v,u,e = log_e(s)
    d15 = all(e[i]==e[i+15] for i in range(len(e)-20))
    print(f"{nom}: v={v} u={u} davr15={'HA ✓' if d15 else 'yoq'}  e1..15={e[:15]}")
    return d15, e
A1=theta(7,8); B1=[0]*(N+1); t=theta(2,13)
for i,x in enumerate(t):
    if i+1<=N: B1[i+1]=x
A2=theta(4,11); B2=[0]*(N+1); t=theta(1,14)
for i,x in enumerate(t):
    if i+1<=N: B2[i+1]=x
for (nom,A,B,ish) in [("(7,8) A^2+AB+B^2",A1,B1,+1), ("(7,8) A^2-AB+B^2",A1,B1,-1),
                      ("(4,11) A^2-AB+B^2",A2,B2,-1), ("(4,11) A^2+AB+B^2",A2,B2,+1)]:
    S=[x+ish*z+y for x,z,y in zip(kop(A,A),kop(A,B),kop(B,B))]
    sinov(nom,S)
