# Знак функционального уравнения w(E1/k_s) для многих сечений; ищем закономерность.
import functools, sys
print = functools.partial(print, flush=True)
def sqf(n):
    n = ZZ(n); return sign(n)*prod(p^(e % 2) for p, e in n.abs().factor())
# примитивные сечения b^2+h^2=2n^2, b>h>0, gcd=1, n <= NMAX
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
secs = []
for n in range(5, NMAX + 1, 2):
    for h in range(1, n, 2):
        b2 = 2*n*n - h*h
        if b2 <= h*h: continue
        b = ZZ(b2).isqrt()
        if b*b == b2 and gcd(b, h) == 1 and gcd(b, n) == 1 and b > h:
            secs.append((b, h, n))
print("сечений:", len(secs))
R = PolynomialRing(QQ, 'z'); z = R.gen(); Fr = R.fraction_field()
T = PolynomialRing(QQ, 'T').gen()
res = []
for (b, h, n) in secs:
    try:
        A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2)
        f = T*(T^2 - 1)*(A*T - C)*(C*T - A)
        Fz = R(Fr((1 - z)^6) * Fr(f((1 + z)/(1 - z))))
        c = Fz.leading_coefficient(); q4 = R(Fz/(c*z)); Ap = q4[2]
        beta = QQ(C - A)/(C + A); dk = sqf((C - A)*(C + A))
        k = QuadraticField(dk, 'r'); sb = k(beta).sqrt()
        Ru = PolynomialRing(k, 'u'); u = Ru.gen()
        D = sqf(ZZ(f(QQ(n + h)/(n + b)).numerator() * f(QQ(n + h)/(n + b)).denominator()))
        cub = Ru(D*c*(u^2 + Ap - 2*beta)*(u + 2*sb))
        a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
        E = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2]).global_minimal_model(semi_global=True)
        nf = pari.nfinit(pari(f'y^2 - {dk}'))
        Ep = pari.ellinit([pari(str(x.polynomial('y'))) for x in E.ainvs()], nf)
        w = ZZ(pari.ellrootno(Ep))
        s_par = QQ(n - h)/(n + b)     # параметр сечения
        res.append((b, h, n, dk, w, s_par))
        print(f"({b},{h},{n}) k=Q(sqrt{dk}) w={w:+d} s={s_par}")
    except Exception as ex:
        print(f"({b},{h},{n}) ошибка {type(ex).__name__}: {str(ex)[:50]}")
pm = sum(1 for t in res if t[4] == -1)
print(f"\nИТОГО: w=-1 у {pm} из {len(res)} ({100.0*pm/max(len(res),1):.0f}%)")
