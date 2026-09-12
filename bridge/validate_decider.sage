# -*- coding: utf-8 -*-
# Валидация решающей процедуры disk_solvable (из local_H_16_5.sage)
# против независимого перебора по Z/p^k с тестом Qp.is_square().
R = PolynomialRing(ZZ, 'u'); u = R.gen()


def unit_square_class_ok(hr, e, p):
    if e % 2 == 1:
        return False
    if p == 2:
        return (hr % 8) == 1
    return kronecker(hr, p) == 1


def prim_part(g, p):
    k = min(ZZ(c).valuation(p) for c in g.coefficients() if c != 0)
    return (ZZ(k), g if k == 0 else g.map_coefficients(lambda c: ZZ(c)//p**k))


def disk_solvable(h, e, p, depth, maxdepth):
    step = 8 if p == 2 else p
    und = False
    for r in range(step):
        hr = ZZ(h(r))
        if hr % p != 0:
            if unit_square_class_ok(hr, e, p):
                return True
            continue
        sub = h(r + step*h.parent().gen())
        if sub.is_zero():
            return True
        k, h2 = prim_part(sub, p)
        e2 = (e + k) % 2
        if depth >= maxdepth:
            und = True
            continue
        res = disk_solvable(h2, e2, p, depth+1, maxdepth)
        if res is True:
            return True
        if res is None:
            und = True
    return None if und else False


def decide_Zp(f, p, maxdepth=25):
    """Решает существование t в Z_p с f(t) в (Q_p)^2 (0 разрешён)."""
    k, h = prim_part(f, p)
    return disk_solvable(h, k % 2, p, 0, maxdepth)


def decide(f, p, maxdepth=25):
    """ПОЛНОЕ решение для P^1(Q_p): карта t в Z_p плюс карта в бесконечности
       t = 1/z, z в p*Z_p (обратный многочлен со сдвигом z = p*w)."""
    K = Qp(p, 60)
    if len(f.change_ring(K).roots()) > 0:
        return True                      # точка Вейерштрасса где-то в Q_p
    r1 = decide_Zp(f, p, maxdepth)
    if r1 is True:
        return True
    fr = f.reverse()
    w = fr.parent().gen()
    r2 = decide_Zp(fr(p*w), p, maxdepth)
    if r2 is True:
        return True
    if r1 is None or r2 is None:
        return None
    return False


def brute_Zp(f, p, k):
    """Независимый перебор по ВСЕЙ P^1(Q_p): карта t в Z_p (целые t < p^k)
       плюс карта в бесконечности t = 1/z, z в p*Z_p (обратный многочлен),
       плюс корни f в Q_p (точки Вейерштрасса вне Z_p)."""
    K = Qp(p, 40)
    for t0 in range(p**k):
        if K(f(t0)).is_square():
            return True, ('t=%d' % t0)
    fr = f.reverse()
    for z0 in range(0, p**k, p):
        if K(fr(z0)).is_square():
            return True, ('inf z=%d' % z0)
    if len(f.change_ring(Qp(p, 60)).roots()) > 0:
        return True, 'Qp-корень'
    return False, None


set_random_seed(20260912)
agree = 0
disagree = 0
nfalse = 0
nundec = 0
for trial in range(200):
    p = choice([2, 3, 5, 7, 11, 13])
    while True:
        f = R([ZZ.random_element(-30, 30) for _ in range(7)])
        if f.degree() == 6 and f.is_squarefree():
            break
    d = decide(f, p)
    if d is None:
        nundec += 1
        continue
    k = 12 if p == 2 else (9 if p == 3 else (7 if p == 5 else 5))
    bt, w = brute_Zp(f, p, k)
    if d == bt:
        agree += 1
    else:
        disagree += 1
        print("РАСХОЖДЕНИЕ p=%d f=%s decide=%s brute=%s %s" % (p, f, d, bt, w))
    if d is False:
        nfalse += 1
print("ВАЛИДАЦИЯ: согласий %d, расхождений %d, ответов False %d, не решено %d"
      % (agree, disagree, nfalse, nundec))
