load('/home/kep/magicKube/family/independent_chabauty/chab.sage')
import sys
b, h, n, D, u0, p = sys.argv[1:7]
prec = int(sys.argv[7]) if len(sys.argv) > 7 else 60
M = int(sys.argv[8]) if len(sys.argv) > 8 else 24
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D)); E = S.E
G = S.point_from_u(QQ(u0))
ch = Chab(S, G, ZZ(p), prec, M).setup()
Qp_ = ch.Qp
pp = ZZ(p)
print(f"E = {E.ainvs()}; N = {ch.N}, orders {ch.orders}, v = {ch.v}, scale = {S.scale}")
print(f"torsion x: {[T[0] for T in ch.tors if not T.is_zero()]}")

def ratrec(a, m, B=None):
    """рациональная реконструкция a mod m."""
    a = ZZ(a) % m
    if B is None: B = isqrt(m//2)
    r0, r1, s0, s1 = m, a, 0, 1
    while r1 > B:
        q = r0//r1
        r0, r1 = r1, r0 - q*r1
        s0, s1 = s1, s0 - q*s1
    if s1 == 0 or gcd(r1, s1) != 1: return None
    return QQ(r1)/QQ(s1) * (1 if s1 > 0 else 1)

def xbranch(Q, i):
    sgn = (1, -1)[i]
    Xs = ch.xseries(Q)
    Rs = PowerSeriesRing(Qp_, 's', default_prec=M); s = Rs.gen()
    zs = Rs(sum(Qp_(ch.emb(c, sgn))*(ch.l[i]*s)^j for j, c in enumerate(ch.ex.list()))).add_bigoh(M)
    f = Rs(0); zp = Rs(1)
    for j in range(M):
        f += ch.emb(Xs[j], sgn)*zp
        zp = (zp*zs).add_bigoh(M)
    return [f[j] for j in range(M)]

def ev(co, s):
    return sum(co[j]*s^j for j in range(len(co)))

tors = ch.tors
for T in tors:
    Tn = 'O' if T.is_zero() else str(T[0])
    for m0 in range(ch.N):
        Q = m0*G + T
        co, gamma, branches, phico = ch.theta(Q)
        n0, mu, tailok, vals = strassmann(co, ch.v, ch.p, M)
        if n0 == 0: continue
        print(f"--- T={Tn} m0={m0}: Strassmann {n0}, mu={mu}, vals={vals[:10]}, gamma={gamma}")
        Rp = PolynomialRing(Qp_, 'S'); Sv = Rp.gen()
        f = sum(co[j]*Sv^j for j in range(M))
        df = f.derivative()
        seen = []
        for a in range(pp):
            s = Qp_(a)
            ok = True
            for _ in range(60):
                fv, dv = ev(co, s), ev([j*co[j] for j in range(1, M)], s)
                if dv == 0 or dv.valuation() > prec//4: ok = False; break
                s = s - fv/dv
            if ok and ev(co, s).valuation() > prec//2:
                if all((s - t).valuation() < prec//4 for t in seen):
                    seen.append(s)
        print(f"    простых нулей найдено: {len(seen)}")
        for s in seen:
            lift = ZZ(s.lift()) if s.valuation() >= 0 else None
            rs = ratrec(lift, pp^(prec - 8)) if lift is not None else None
            xs = []
            for i in (0, 1):
                try:
                    xs.append(ev(xbranch(Q, i), s))
                except Exception as ex:
                    xs.append(None)
            info = ""
            if xs[0] is not None and xs[0].valuation() >= 0:
                xl = ZZ(xs[0].lift()); rx = ratrec(xl, pp^(prec - 8))
                info = f" x={xs[0] + O(pp^10)} ratrec(x)={rx} u={None if rx is None else rx/S.scale}"
            print(f"    s = {s + O(pp^10)}  ratrec(s) = {rs}{info}")
        if ch.N % 2 == 0 and m0 == ch.N//2:
            hv = ev(co, Qp_(-1)/2)
            print(f"    theta(-1/2): val = {hv.valuation()}")
