\\ ildiz_gp.gp — PARI/GP da NOLDAN (Python kodidan mustaqil) tekshiruv
default(parisize,"1G"); default(nbthreads,1);
N = 1500;
\\ f(-q^a,-q^b) yig'indidan, q — darajali qator
th(a,b) = { my(s = O(q^N)); for(n = -ceil(sqrt(2*N/(a+b)))-3, ceil(sqrt(2*N/(a+b)))+3, my(e = a*n*(n+1)/2 + b*n*(n-1)/2); if(e >= 0 && e < N, s += (-1)^n*q^e)); s };
\\ Jacobi ko'paytmasidan
jp(a,b) = { my(M = a+b, s = 1 + O(q^N)); forstep(k = a, N, M, s *= (1 - q^k)); forstep(k = b, N, M, s *= (1 - q^k)); forstep(k = M, N, M, s *= (1 - q^k)); s };
\\ Borwein a(q^k)
ab(k) = { my(s = O(q^N), M = ceil(sqrt(4*N/(3*k)))+3); for(m = -M, M, for(n = -M, M, my(e = k*(m^2+m*n+n^2)); if(e < N, s += q^e))); s };
ok = 0; jami = 0;
chk(nom, L, R) = { jami++; if(L == R, ok++; print("  MOS  ", nom), print("  FARQ ", nom, "  birinchi farq: ", valuation(L-R,q))); };
print("== G1: yig'indi == Jacobi ko'paytma");
foreach([[1,4],[2,3],[3,12],[6,9],[7,8],[2,13],[4,11],[1,14],[5,10]], v, chk(Str("f(-q^",v[1],",-q^",v[2],")"), th(v[1],v[2]), jp(v[1],v[2])));
A1 = th(7,8); B1 = q*th(2,13); A2 = th(4,11); B2 = q*th(1,14);
f14 = th(1,4); f23 = th(2,3); f312 = th(3,12); f69 = th(6,9); a5 = ab(5);
print("== G2: Theorem 1, 2 (", N, " had)");
chk("Theorem 1:  A^3-B^3 = q^2 f312^3 + f23 a(q^5)", A1^3 - B1^3, q^2*f312^3 + f23*a5);
chk("Theorem 2:  q(A^3+B^3) = f69^3 - f14 a(q^5)", q*(A2^3 + B2^3), f69^3 - f14*a5);
chk("by-product 8.2.1: 3AB(A-B) f312 = f23 (a(q^5) f312 - f14^3)", 3*A1*B1*(A1-B1)*f312, f23*(a5*f312 - f14^3));
chk("by-product 8.2.2: q 3AB(A+B) f69 = f14 (a(q^5) f69 - f23^3)", q*3*A2*B2*(A2+B2)*f69, f14*(a5*f69 - f23^3));
print("== G3: Borwein a^3 = b^3 + c^3 (bo'lishsiz)");
E1 = prod(k=1, N, 1 - q^k + O(q^N)); E3 = prod(k=1, N\3, 1 - q^(3*k) + O(q^N));
chk("a^3 E1^3 E3^3 = E1^12 + 27 q E3^12", ab(1)^3*E1^3*E3^3, E1^12 + 27*q*E3^12);
print("== G4: manfiy nazorat (FARQ kutiladi)");
jami++; if(A1^3 - B1^3 != q^2*f312^3 + f23*(a5 + q^10), ok++; print("  MOS   nazorat FARQ berdi (sezgir)"), print("  FARQ  nazorat sezmadi!"));
print("\nGP-ILDIZ: ", ok, "/", jami, " MOS");
quit
