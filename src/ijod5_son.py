# N1, N2 ni mpmath 50 xonada, a(q) ni Entry 8.2.3 dagi {b^3+c^3}^(1/3) shaklidan EMAS, to'g'ridan qatordan hisoblab
from mpmath import mp, mpf
mp.dps=50
def f(a,b):
    s=mpf(0); n=0
    while True:
        t1=a**(n*(n+1)//2)*b**(n*(n-1)//2); t2=(a**((-n)*(-n+1)//2)*b**((-n)*(-n-1)//2)) if n else 0
        s+=t1+t2
        if abs(t1)+abs(t2)<mpf(10)**-60 and n>3: break
        n+=1
    return s
def a_(x):
    s=mpf(0); M=1
    while x**(M*M*3//4) > mpf(10)**-60: M+=1
    for m in range(-M,M+1):
        for n in range(-M,M+1): s+=x**(m*m+m*n+n*n)
    return s
w=mpf(0)
for q in [mpf('0.1'),mpf('0.3'),mpf('0.5'),mpf('0.7')]:
    A1,B1=f(-q**7,-q**8),q*f(-q**2,-q**13); A2,B2=f(-q**4,-q**11),q*f(-q,-q**14)
    e1=abs((A1**3-B1**3)-(q**2*f(-q**3,-q**12)**3+f(-q**2,-q**3)*a_(q**5)))
    e2=abs((A2**3+B2**3)-(f(-q**6,-q**9)**3-f(-q,-q**4)*a_(q**5))/q)
    w=max(w,e1,e2); print('q=',q,' N1 farq',mp.nstr(e1,3),' N2 farq',mp.nstr(e2,3))
print('JAMI:', mp.nstr(w,3), '->', 'MOS' if w<mpf(10)**-40 else 'FARQ')
