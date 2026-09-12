# Независимая перепроверка local_H_11_4.sage другим кодом и другой логикой.
#  (1) sanity-тесты для Rational.is_padic_square;
#  (2) тупой перебор целых t в [-P, P] + t = 1/u, точный тест f(t) in (Q_p^*)^2;
#  (3) для p = 1 mod 4 — явный корень секстики в Z_p (точка Вейерштрасса Y = 0),
#      найденный через Hensel-подъём sqrt(-1), с точной верификацией;
#  (4) счёт точек редукции над F_p для хороших простых (граница Вейля).
RZ = PolynomialRing(ZZ, 't'); T = RZ.gen()
f = RZ(530464*T^6 + 4612242*T^4 + 4612242*T^2 + 530464)
assert f == RZ(4*(16*T^2+121)*QQ(137)/2*(1+T^2)*(121*T^2+16)*1), "модель"
print("f =", factor(f))

print("\n(1) sanity is_padic_square:")
for (x, p, exp) in [(9,7,True), (7,7,False), (2,7,True), (3,7,False),
                    (17,2,True), (3,2,False), (2,2,False), (QQ(1)/4,2,True),
                    (137,2,True), (8,2,False), (4,2,True), (49*5,5,False), (25,5,True)]:
    got = QQ(x).is_padic_square(p)
    print("   %-8s p=%-4d ожидалось %-6s получено %-6s %s" % (x, p, exp, got, "OK" if got==exp else "!!!"))

badp = [2,3,5,7,11,137]
allp = sorted(set(badp) | set(primes(2,60)))

print("\n(2) тупой перебор: минимальное |t| (целое) или u=1/t с f — квадратом в Q_p")
res2 = {}
for p in allp:
    found = None
    for k in range(0, 20000):
        for tv in ([0] if k == 0 else [k, -k]):
            v = QQ(f(tv))
            if v == 0 or v.is_padic_square(p):
                found = ("t=%s" % tv, v); break
        if found: break
    if not found:
        for k in range(1, 20000):
            for uv in [k, -k]:
                v = QQ(f(QQ(1)/uv))
                if v == 0 or v.is_padic_square(p):
                    found = ("t=1/%s" % uv, v); break
            if found: break
    res2[p] = found
    print("   p=%-5d %s" % (p, "%s,  f = %s" % found if found else "НЕ НАЙДЕНО в диапазоне"))

print("\n(3) p = 1 mod 4: корень t^2+1 в Z_p (Hensel), точка (t,0) на H")
res3 = {}
for p in [q for q in allp if q % 4 == 1]:
    # подъём sqrt(-1) mod p^40
    r = ZZ(Mod(-1, p).sqrt().lift())
    N = 1
    while N < 40:
        pk = p**N
        # Ньютон: r <- r - (r^2+1)/(2r) mod p^(2N)
        M = p**(2*N)
        r = ZZ(Mod(r - (r*r+1)*ZZ(Mod(2*r, M)**(-1)), M))
        N = 2*N
    M = p**N
    ok = (r*r + 1) % M == 0
    fv = ZZ(f(r)) % M
    res3[p] = ok and fv % p**min(N, 30) == 0
    print("   p=%-5d  r^2+1 = 0 mod p^%d : %s;  f(r) = 0 mod p^%d : %s  => точка (t,Y)=(r,0)"
          % (p, N, ok, min(N,30), fv % p**min(N,30) == 0))

print("\n(4) хорошие простые: число точек гладкой редукции H над F_p и граница Вейля")
for p in [q for q in primes(2, 200) if q not in badp]:
    Fp = GF(p)
    cnt = 0
    for x in Fp:
        v = Fp(f(ZZ(x)))
        cnt += (1 if v == 0 else (2 if v.is_square() else 0))
    lcp = Fp(f.leading_coefficient())
    cnt += (2 if lcp.is_square() else 0)       # точки на бесконечности (lc != 0 при хорошей редукции)
    weil = p + 1 - 4*RR(p).sqrt()
    if p < 70 or cnt == 0:
        print("   p=%-5d #H(F_p) = %-5d  граница Вейля p+1-4sqrt(p) = %.2f  %s"
              % (p, cnt, weil, "" if cnt > 0 else "<-- ПУСТО!"))
    assert cnt > 0, "пустая редукция при p=%d" % p
print("   все хорошие простые p < 200: #H(F_p) > 0 (гладкая точка => подъём по Гензелю)")
print("   граница Вейля p+1-4*sqrt(p) > 0 при p >= 17 => для всех хороших p >= 17 точка есть безусловно")
