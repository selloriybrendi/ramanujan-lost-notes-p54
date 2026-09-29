# IJOD-1: yo'qolgan tomondagi kvadratik omil Q = A^2 + AB + B^2 (8.2.1) va A^2 - AB + B^2 (8.2.2).
# G'oya: A, B — bir xil q^(1/120) siljishli vazn-1/2 theta => Q vazn-1 modulyar shakl.
# Vazn-1 shakllar ko'pincha Eisenstein (Lambert qator) => koeffitsiyentlar KICHIK va "bo'luvchi-yig'indi" tuzilishli.
N = 600
def theta(a,b):
    s=[0]*(N+1); n=0
    while True:
        did=False
        for m in ((n,-n) if n else (0,)):
            e=a*m*(m+1)//2+b*m*(m-1)//2
            if 0<=e<=N: s[e]+=(-1)**(m%2); did=True
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
def sil(s,k): return [0]*k+s[:N+1-k]
A1,B1=theta(7,8),sil(theta(2,13),1)
A2,B2=theta(4,11),sil(theta(1,14),1)
Q1=[a+m+b for a,m,b in zip(kop(A1,A1),kop(A1,B1),kop(B1,B1))]
Q2=[a-m+b for a,m,b in zip(kop(A2,A2),kop(A2,B2),kop(B2,B2))]
print("Q1 (8.2.1, A2+AB+B2):", Q1[:61])
print("Q2 (8.2.2, A2-AB+B2):", Q2[:61])
print("max|Q1| n<=600:", max(abs(x) for x in Q1), " max|Q2|:", max(abs(x) for x in Q2))
