# Разбор единственного расхождения из validate_decider.sage
R = PolynomialRing(ZZ, 'u'); u = R.gen()
f = 20*u^6 - 21*u^5 + 22*u^4 - 23*u^3 + 3*u^2 + 22*u - 21
p = 2
K = Qp(2, 60)

print("f =", f)
print("f(r) mod 8 для r=0..7:", [(r, ZZ(f(r)) % 8, ZZ(f(r)).valuation(2)) for r in range(8)])

# развернём рекурсию вручную, ища сертификат
def prim_part(g, p):
    k = min(ZZ(c).valuation(p) for c in g.coefficients() if c != 0)
    return (ZZ(k), g if k == 0 else g.map_coefficients(lambda c: ZZ(c)//p**k))

def search(h, e, off, scale, depth, maxdepth=25):
    """h -- примитивный; f(off + scale*u) = 2^(что-то с чётностью e) * h(u)"""
    for r in range(8):
        hr = ZZ(h(r))
        if hr % 2 != 0:
            if e % 2 == 0 and hr % 8 == 1:
                return (off + scale*r, scale*8, depth)
            continue
        sub = h(r + 8*h.parent().gen())
        if sub.is_zero():
            return (off + scale*r, scale*8, depth)
        k, h2 = prim_part(sub, 2)
        if depth >= maxdepth:
            continue
        res = search(h2, (e+k) % 2, off + scale*r, scale*8, depth+1, maxdepth)
        if res:
            return res
    return None

k0, h0 = prim_part(f, 2)
cert = search(h0, k0 % 2, 0, 1, 0)
print("сертификат (центр диска, шаг, глубина):", cert)
if cert:
    t0, step, d = cert
    print("v_2(шага) =", ZZ(step).valuation(2))
    for j in range(6):
        tt = t0 + j*step
        v = ZZ(f(tt))
        print("   t=%-14d f(t)=%-24d v_2=%-3d единица mod 8 = %-2d  квадрат в Q_2: %s"
              % (tt, v, ZZ(v).valuation(2), ZZ(v // 2**ZZ(v).valuation(2)) % 8, K(v).is_square()))
    print("=> точка существует; грубый перебор по целым t < 2^12 её не видел,")
    print("   потому что диск имеет уровень 2^%d." % ZZ(step).valuation(2))
