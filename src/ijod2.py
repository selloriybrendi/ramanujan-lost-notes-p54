# IJOD-2: T(x) = sum Q[m] x^(60m+1) — Hecke xos-shaklmi? multiplikativlik va tub a_p lar
from math import gcd
import ijod1 as I
def T_of(Q):
    return {60*m+1: c for m, c in enumerate(Q)}
def tub(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
for nom,Q in (("Q1 (8.2.1)",I.Q1),("Q2 (8.2.2)",I.Q2)):
    T=T_of(Q); keys=sorted(T)
    ok=bad=0; misol=[]
    for m in keys:
        for n in keys:
            if m<n and gcd(m,n)==1 and m*n in T and m>1:
                if T[m]*T[n]==T[m*n]: ok+=1
                else: bad+=1; misol.append((m,n,T[m],T[n],T[m*n]))
    print(f"{nom}: multiplikativ sinov: mos {ok}, buzilgan {bad}", misol[:3])
    ps=[(p,T[p]) for p in keys if tub(p)][:25]
    print("   a_p (p tub, p=60m+1):", ps)
