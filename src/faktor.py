# A^3-B^3 = (A-B)(A^2+AB+B^2). A-B ma'lum (8.2.10). A^2+AB+B^2 ham mahsulotmi?
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
    num = [i*c[i] for i in range(len(c))]; M = len(c)-1; L = [0]*(M+1)
    for m in range(M+1):
        t = num[m]
        for k in range(1, m+1): t -= c[k]*L[m-k] if k < len(c) else 0
        L[m] = t//c[0]
    e = [0]*(M+1)
    for m in range(1, M+1):
        t = -L[m]
        for d in range(1, m):
            if m % d == 0: t -= d*e[d]
        e[m] = t//m
    return e[1:]
A = theta(4,11); Bq = theta(1,14)
B = [0]*(N+1)
for i,x in enumerate(Bq):
    if i+1 <= N: B[i+1] = x
Q = [a*a for a in [0]]  # placeholder
A2 = kop(A,A); B2 = kop(B,B); AB = kop(A,B)
S = [x+y+z for x,y,z in zip(A2,AB,B2)]
e = log_e(S)
print("A^2+AB+B^2 e_1..30:", e[:30])
print("davr15:", all(e[i]==e[i+15] for i in range(len(e)-20)))
print("A^2+AB+B^2 koef 0..24:", S[:25])
# A^3+B^3 ham (yo'lda):
A3B3p = [x+y for x,y in zip(kop(A2,A), kop(B2,B))]
ep = log_e(A3B3p)
print("A^3+B^3 e_1..30:", ep[:30], "davr15:", all(ep[i]==ep[i+15] for i in range(len(ep)-20)), "davr30:", all(ep[i]==ep[i+30] for i in range(len(ep)-35)))
# OEIS uchun nishon koef:
A3B3 = [x-y for x,y in zip(kop(A2,A), kop(B2,B))]
print("A^3-B^3 koef 0..20:", A3B3[:21])
