# Fable, 25.09.2026. Независимая проверка certificate_<p>_<q>.json: отображения, формулы сложения, тождества Безу,
# нормы, точный экспоненциальный сертификат. Строится заново из наборов клеток; читаются только коэффициенты.
import sys, json
sys.path.insert(0, '/home/kep/magicKube/rigorous_fable')
from fmodels import *

def dec(V, data):
    vs = V.gens(); out = V(0)
    for e, c in data:
        c = QQ(c); assert c.denominator() == 1
        out += c*prod(v**int(k) for v, k in zip(vs, e))
    return out

def check_P1(c, phi_expected, tvar):
    R = tvar.parent(); t = tvar; K = R.fraction_field()
    F = R(list(map(QQ, c['F']))); G = R(list(map(QQ, c['G']))); d = c['degree']; r = ZZ(c['R'])
    assert r != 0 and max(F.degree(), G.degree()) == d and F.gcd(G) == 1
    assert all(v.denominator() == 1 for v in F.list() + G.list())
    assert K(F/G) == phi_expected
    norms = []
    for w in c['witnesses']:
        A = R(list(map(QQ, w['A']))); B = R(list(map(QQ, w['B'])))
        assert all(v.denominator() == 1 for v in A.list() + B.list()) and A.degree() <= d - 1 and B.degree() <= d - 1
        assert A*F + B*G == r*t**w['j']; norms.append(sum(abs(v) for v in A.list() + B.list()))
    assert sorted(w['j'] for w in c['witnesses']) == [0, 2*d - 1]
    assert ZZ(c['K']) == max(norms) and ZZ(c['L']) == max(sum(abs(v) for v in F.list()), sum(abs(v) for v in G.list()))

def verify(fn):
    Dt = json.load(open(fn)); p, q = map(ZZ, Dt['slope'].split('/'))
    assert [QQ(c) for c in Dt['cells']] == cells(p, q)
    Ts = [[QQ(c) for c in T] for T in Dt['Ts']]; hub = Dt['hub']
    Et = EllipticCurve(QQ, list(map(QQ, Dt['E']))); assert Et.a1() == 0 and Et.a3() == 0
    facs = [Fac(T) for T in Ts]; ims = [Image(f, Et) for f in facs]   # конструкторы проверяют тождества кривых
    for im, rec in zip(ims, Dt['images']):
        assert im.M == matrix(ZZ, [[ZZ(c) for c in row] for row in rec['M']]) and im.mobius_norm() == ZZ(rec['mobius_norm'])
    Rt = PolynomialRing(QQ, 't'); t = Rt.gen(); Ft = Rt.fraction_field()
    b2, b4, b6, b8 = Et.b_invariants()
    check_P1(Dt['duplication'], Ft((t**4 - b4*t*t - 2*b6*t - b8)/(4*t**3 + b2*t*t + 2*b4*t + b6)), t)
    Bv = ZZ(Dt['B']); assert Bv == max(ZZ(Dt['duplication']['K']), ZZ(Dt['duplication']['L']))
    for c in Dt['pairs']:
        i, j = c['i'], c['j']; assert i == hub
        A, B, D, I = pair_AB(ims[i], ims[j]); assert [str(x) for x in D] == c['D']
        assert exact_pair_check(ims[i], ims[j], A, B, D, QQ(1)/(7*max(abs(x) for x in cells(p, q))))
        assert exact_pair_check(ims[i], ims[j], A, B, D, QQ(-1)/(13*max(abs(x) for x in cells(p, q))))
        if c['kind'] == 'conic':
            c1, c2 = map(QQ, c['D']); tt = Ft(t)
            if c1 == -c2: z = 2*tt/(c2*(1 + tt**2)); y = (1 - tt**2)/(1 + tt**2)
            else: z = (c1 + c2 - 2*tt)/(tt**2 - c1*c2); y = 1 + tt*z
            assert y**2 == (1 + c1*z)*(1 + c2*z)
            Az = rational_function_eval(A, z); Bz = rational_function_eval(B, z); up, um = Az + Bz*y, Az - Bz*y
            if c['compressed']:
                assert c['compressed'] == 'w=t^2'
                def comp(phi):
                    F, G = phi.numerator(), phi.denominator()
                    assert all(F[k] == 0 and G[k] == 0 for k in range(1, max(F.degree(), G.degree()) + 1, 2))
                    return Ft(Rt([F[2*k] for k in range(F.degree()//2 + 1)])/Rt([G[2*k] for k in range(G.degree()//2 + 1)]))
                up, um = comp(up), comp(um)
            check_P1(c['plus'], up, t); check_P1(c['minus'], um, t)
        else:
            aux = Aux(A, B, D); assert [str(a) for a in aux.E3.ainvs()] == c['E3'] and c['basis'] == '1,x,y,x^2'
            basis = aux.basis_3O_plus_S(None)
            V = PolynomialRing(QQ, ['v0', 'v1', 'v2', 'v3']); vs = V.gens()
            quads = [dec(V, qd) for qd in c['quadrics']]
            for qd in quads:   # квадрики обращаются в нуль на кривой
                assert qd.is_homogeneous() and qd.degree() == 2
                assert sum(coef*prod(phi**int(k) for phi, k in zip(basis, e)) for e, coef in qd.dict().items()) == 0
            for name, u in (('plus', aux.up), ('minus', aux.um)):
                rec = c[name]; n = list(map(QQ, rec['n'])); h = list(map(QQ, rec['h']))
                assert sum(a*phi for a, phi in zip(n, basis)) == u*sum(a*phi for a, phi in zip(h, basis))
                F = sum(a*v for a, v in zip(n, vs)); G = sum(a*v for a, v in zip(h, vs)); N = rec['N']; r = ZZ(rec['R']); assert r != 0
                norms = []
                for v, w in zip(vs, rec['witnesses']):
                    pols = [dec(V, pd) for pd in w['polynomials']]; assert len(pols) == 2 + len(quads)
                    assert all(pp == 0 or (pp.is_homogeneous() and pp.degree() == N - 1) for pp in pols[:2])
                    assert all(pp == 0 or (pp.is_homogeneous() and pp.degree() == N - 2) for pp in pols[2:])
                    assert sum(pp*g for pp, g in zip(pols, [F, G] + quads)) == r*v**N
                    norms.append(sum(abs(cc) for pp in pols[:2] for cc in pp.coefficients()))
                assert len(rec['witnesses']) == 4 and ZZ(rec['K']) == max(norms)
                assert ZZ(rec['L']) == max(sum(abs(cc) for cc in F.coefficients()), sum(abs(cc) for cc in G.coefficients()))
        Av = ZZ(c['A']); assert Av == max(ZZ(c['plus']['L'])*ZZ(c['minus']['K']), ZZ(c['minus']['L'])*ZZ(c['plus']['K']))
        C = ZZ(c['integer_bound']); partial = QQ(1); term = QQ(1)
        for k in range(1, c['exp_terms'] + 1):
            term *= QQ(12*C)/k; partial += term
        assert partial > Av**3*Bv**2
        print(fn, 'pair', c['i'], c['j'], c['kind'], 'verified; |<Q_i,Q_j>| <', C, flush=True)
    print(fn, 'all identities verified; |hhat - h(x)| <= log(B)/3, B =', Bv, flush=True)

if __name__ == '__main__':
    for fn in sys.argv[1:]: verify(fn)
