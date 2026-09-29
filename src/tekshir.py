# Gipoteza: A^3 - B^3 = (q^5;q^5)_inf^3 * f(-q^3,-q^12) / f(-q^6,-q^9)
# Ekvivalent Euler-ko'rsatkichlar (mod 15): n=3,12:+1 · n=6,9:-1 · n=5,10,0:+3
N = 400
def poch(exps):  # prod (1-q^n)^{e(n mod 15)}
    s = [0]*(N+1); s[0] = 1
    for n in range(1, N+1):
        e = exps[n % 15]
        if e == 0: continue
        # (1-q^n)^e ga ko'paytirish
        for _ in range(abs(e)):
            if e > 0:
                for i in range(N, n-1, -1): s[i] -= s[i-n]
            else:
                for i in range(n, N+1): s[i] += s[i-n]
    return s
exps = {0:3, 1:0, 2:0, 3:1, 4:0, 5:3, 6:-1, 7:0, 8:0, 9:-1, 10:3, 11:0, 12:1, 13:0, 14:0}
G = poch(exps)

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
A = theta(4, 11); Bq = theta(1, 14)
B = [0]*(N+1)
for i, x in enumerate(Bq):
    if i+1 <= N: B[i+1] = x
A3B3 = [x-y for x, y in zip(kop(kop(A,A),A), kop(kop(B,B),B))]
mos = A3B3 == G
print("A^3-B^3 == (q^5;q^5)^3 * f(-q^3,-q^12)/f(-q^6,-q^9) [400 hadgacha]:", "MOS ✓" if mos else "FARQ ✗")
if not mos:
    d = next(i for i in range(N+1) if A3B3[i] != G[i]); print("birinchi farq q^", d, A3B3[d], G[d])
# Qo'shimcha: A^3+B^3 ham davriymi?
import subprocess
