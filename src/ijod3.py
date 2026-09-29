# IJOD-3: n=60m+1 ning binar formalar bilan tasvirlari va Q-koeffitsiyent yonma-yon
import ijod1 as I
T1={60*m+1:c for m,c in enumerate(I.Q1)}; T2={60*m+1:c for m,c in enumerate(I.Q2)}
def rep(n,a,b,c):
    out=[]
    X=int((4*a*n)**0.5)+3
    for x in range(-X,X+1):
        for y in range(0,X+1):
            if a*x*x+b*x*y+c*y*y==n: out.append((x,y))
    return out
for n in sorted(T1)[:14]:
    print(n, "Q1:",T1[n],"Q2:",T2[n], "| x2+15y2:",[(x,y) for x,y in rep(n,1,0,15) if x>=0], " 3x2+5y2:",[(x,y) for x,y in rep(n,3,0,5) if x>=0],
          " x2+xy+4y2(D-15):", len(rep(n,1,1,4)), " 2x2+xy+2y2:", len(rep(n,2,1,2)))
