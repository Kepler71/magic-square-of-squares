# Fable, 25.09.2026. Сертификаты спариваний: тождества Безу для u± (коника: P^1; эллиптическая: вложение [1:x:y:x^2] в P^3)
# и для формулы удвоения на E. Выход: certificate_<p>_<q>.json. Ничего из кода Codex не импортируется.
import sys, json, time
sys.path.insert(0, '/home/kep/magicKube/rigorous_fable')
from fmodels import *

def enc(p):
    return [[list(map(int, e)), str(c)] for e, c in p.dict().items()]

def l1(p):
    return ZZ(sum(abs(c) for c in p.coefficients())) if hasattr(p, 'coefficients') else ZZ(sum(abs(c) for c in p))

def cert_P1(phi, tvar):
    """phi ∈ Q(t); F,G ∈ Z[t] взаимно простые; d = max deg; тождества A F + B G = R t^j, j ∈ {0, 2d−1}, deg A,B ≤ d−1."""
    R = tvar.parent(); t = tvar
    F, G = R(phi.numerator()), R(phi.denominator())
    den = lcm([c.denominator() for c in F.list() + G.list()]); F, G = R(F*den), R(G*den)
    cont = gcd([ZZ(c) for c in F.list() + G.list()]); F, G = R(F/cont), R(G/cont)
    assert F.gcd(G) == 1
    d = max(F.degree(), G.degree())
    def vec(P): l = P.list(); return l + [QQ(0)]*(2*d - len(l))
    S = matrix(QQ, [vec(t**i*H) for H in (F, G) for i in range(d)]).transpose(); assert S.det() != 0
    sols = {}
    for j in (0, 2*d - 1):
        e = vector(QQ, [int(i == j) for i in range(2*d)]); sols[j] = S.solve_right(e)
    Rres = ZZ(lcm([c.denominator() for s in sols.values() for c in s]))
    wits = []
    for j, s in sols.items():
        s = s*Rres; assert all(c.denominator() == 1 for c in s)
        A = R(list(s[:d])); B = R(list(s[d:])); assert A*F + B*G == Rres*t**j
        wits.append({'j': int(j), 'A': [str(c) for c in A.list()], 'B': [str(c) for c in B.list()], 'norm': str(l1(A) + l1(B))})
    return {'degree': int(d), 'F': [str(c) for c in F.list()], 'G': [str(c) for c in G.list()], 'R': str(Rres),
            'K': str(max(ZZ(w['norm']) for w in wits)), 'L': str(max(l1(F), l1(G))), 'witnesses': wits}

def conic_certificate(A, B, D):
    """u±(t) на конике y_D^2=(1+c1 z)(1+c2 z): параметризация через (z,y)=(0,1) (для c1=−c2 — как у Codex: z=2t/(c(1+t^2)))."""
    c1, c2 = sorted(QQ(c) for c in D); Rt = PolynomialRing(QQ, 't'); t = Rt.gen(); Ft = Rt.fraction_field(); t = Ft(t)
    if c1 == -c2:
        c = c2; z = 2*t/(c*(1 + t**2)); y = (1 - t**2)/(1 + t**2); param = 'z=2t/(%s(1+t^2)), y=(1-t^2)/(1+t^2)' % c
    else:
        z = (c1 + c2 - 2*t)/(t**2 - c1*c2); y = 1 + t*z; param = 'z=(c1+c2-2t)/(t^2-c1c2), y=1+tz'
    assert y**2 == (1 + c1*z)*(1 + c2*z)
    Az = rational_function_eval(A, z); Bz = rational_function_eval(B, z)
    up, um = Az + Bz*y, Az - Bz*y
    out = {'kind': 'conic', 'D': [str(c) for c in (c1, c2)], 'param': param, 'up': str(up), 'um': str(um)}
    def even(phi):
        F, G = phi.numerator(), phi.denominator()
        return all(F[i] == 0 and G[i] == 0 for i in range(1, max(F.degree(), G.degree()) + 1, 2))
    if even(up) and even(um):
        def comp(phi):
            F, G = phi.numerator(), phi.denominator()
            return Ft(Rt([F[2*i] for i in range(F.degree()//2 + 1)]) / Rt([G[2*i] for i in range(G.degree()//2 + 1)]))
        up, um = comp(up), comp(um); out['compressed'] = 'w=t^2'
    else: out['compressed'] = None
    cp, cm = cert_P1(up, Rt.gen()), cert_P1(um, Rt.gen()); assert cp['degree'] == cm['degree']
    out['plus'], out['minus'] = cp, cm
    out['A'] = str(max(ZZ(cp['L'])*ZZ(cm['K']), ZZ(cm['L'])*ZZ(cp['K'])))
    return out

def bezout_P3(F, G, quadrics, N=3):
    """тождества A_i F + B_i G + Σ C_ij Q_j = R v_i^N (deg A,B = N−1, deg C = N−2); возвращает R и свидетелей."""
    V = F.parent(); vs = V.gens(); gens = [F, G] + list(quadrics)
    def mons(d): return [prod(v**int(e) for v, e in zip(vs, ee)) for ee in IntegerVectors(d, len(vs))]
    mm = mons(N)
    basis = [(j, m) for j, g in enumerate(gens) for m in mons(N - g.degree())]
    cols = [m*gens[j] for j, m in basis]
    M = matrix(QQ, [[p.monomial_coefficient(m) for p in cols] for m in mm])
    if M.rank() < len(mm): return None
    sols = [M.solve_right(vector(QQ, [(v**N).monomial_coefficient(m) for m in mm])) for v in vs]
    Rres = ZZ(lcm(c.denominator() for s in sols for c in s)); wits = []
    for v, s in zip(vs, sols):
        pols = [V(0)]*len(gens)
        for coef, (j, m) in zip(s, basis): pols[j] += Rres*coef*m
        assert all(c.denominator() == 1 for p in pols for c in p.coefficients())
        assert sum(p*g for p, g in zip(pols, gens)) == Rres*v**N
        wits.append({'polynomials': [enc(p) for p in pols], 'norm': str(l1(pols[0]) + l1(pols[1]))})
    return Rres, wits

def elliptic_certificate(A, B, D):
    aux = Aux(A, B, D); basis = aux.basis_3O_plus_S(None)
    out = {'kind': 'elliptic', 'D': [str(c) for c in aux.D], 'E3': [str(a) for a in aux.E3.ainvs()], 'basis': '1,x,y,x^2'}
    V, quads = vanishing_forms(basis, 2); assert len(quads) == 2
    out['quadrics'] = [enc(q) for q in quads]
    vs = V.gens()
    for name, u in (('plus', aux.up), ('minus', aux.um)):
        sols = solve_linear_forms(u, basis); assert len(sols) == 1
        n, h = sols[0]; sc = lcm([c.denominator() for c in list(n) + list(h)]); n, h = n*sc, h*sc
        g = gcd([ZZ(c) for c in list(n) + list(h)]); n, h = n/g, h/g
        assert sum(c*phi for c, phi in zip(n, basis)) == u*sum(c*phi for c, phi in zip(h, basis))
        F = sum(QQ(c)*v for c, v in zip(n, vs)); G = sum(QQ(c)*v for c, v in zip(h, vs))
        for N in (3, 4):
            res = bezout_P3(F, G, quads, N)
            if res is not None: break
        assert res is not None
        Rres, wits = res
        out[name] = {'n': [str(c) for c in n], 'h': [str(c) for c in h], 'N': int(N), 'R': str(Rres), 'witnesses': wits,
                     'K': str(max(ZZ(w['norm']) for w in wits)), 'L': str(max(l1(F), l1(G)))}
    out['A'] = str(max(ZZ(out['plus']['L'])*ZZ(out['minus']['K']), ZZ(out['minus']['L'])*ZZ(out['plus']['K'])))
    return out

def duplication_certificate(E):
    Rt = PolynomialRing(QQ, 't'); t = Rt.gen(); b2, b4, b6, b8 = E.b_invariants()
    dup = (t**4 - b4*t*t - 2*b6*t - b8)/(4*t**3 + b2*t*t + 2*b4*t + b6)
    c = cert_P1(dup, t); assert c['degree'] == 4
    c['B'] = str(max(ZZ(c['K']), ZZ(c['L']))); return c

def exp_certificate(C, target, maxterms=20000):
    """наименьшее J с Σ_{j≤J} (12C)^j/j! > target (точно над Q); тогда 12C > log target."""
    partial = QQ(1); term = QQ(1)
    for j in range(1, maxterms):
        term *= QQ(12*C)/j; partial += term
        if partial > target: return j
    raise AssertionError('exp certificate failed')

def certify(p, q, Ts, hub, others, tag=None):
    t0 = time.time()
    facs = [Fac(T) for T in Ts]; Et = target_model(facs[hub].E)
    for f in facs: assert f.E.is_isomorphic(Et)
    ims = [Image(f, Et) for f in facs]
    out = {'slope': '%d/%d' % (p, q), 'cells': [str(c) for c in cells(p, q)], 'Ts': [[str(c) for c in T] for T in Ts], 'hub': hub, 'others': others,
           'E': [str(a) for a in Et.ainvs()], 'conductor': str(Et.conductor()),
           'images': [{'T': [str(c) for c in im.fac.T], 'M': [[str(c) for c in row] for row in im.M.rows()], 'mobius_norm': str(im.mobius_norm())} for im in ims]}
    dup = duplication_certificate(Et); out['duplication'] = dup; Bv = ZZ(dup['B'])
    I200 = RealIntervalField(200); pairs = []
    for j in others:
        A, B, D, I = pair_AB(ims[hub], ims[j])
        zt = QQ(1)/(10*max(abs(c) for c in cells(p, q)))
        assert exact_pair_check(ims[hub], ims[j], A, B, D, zt)
        c = conic_certificate(A, B, D) if len(D) == 2 else elliptic_certificate(A, B, D)
        c['i'], c['j'], c['I'] = hub, j, [str(x) for x in I]
        Av = ZZ(c['A']); bnd = (3*I200(Av).log() + 2*I200(Bv).log())/12; C = ZZ(bnd.upper().ceil())
        J = exp_certificate(C, Av**3*Bv**2)
        c['pairing_bound_interval'] = str(bnd); c['integer_bound'] = int(C); c['exp_terms'] = int(J)
        print('pair', Ts[hub], Ts[j], c['kind'], 'L,K +:', c['plus']['L'], c['plus']['K'], '-:', c['minus']['L'], c['minus']['K'], 'bound', bnd, '->', C, 'secs', time.time() - t0, flush=True)
        pairs.append(c)
    out['pairs'] = pairs; out['B'] = str(Bv); out['seconds'] = time.time() - t0
    fn = 'certificate_%s.json' % (tag or '%d_%d' % (p, q)); json.dump(out, open(fn, 'w'), indent=1)
    print('written', fn, 'B', Bv, 'log B/3', I200(Bv).log()/3, flush=True)
    return out

CONFIGS = {
    '126/451': dict(p=126, q=451, Ts=[[-577,-451,-325,-126],[-325,-126,126],[-126,126,325]], hub=1, others=[2, 0]),
    '73/362':  dict(p=73, q=362, Ts=[[-435,-362,-289,-73],[-289,-73,73],[-73,73,289]], hub=1, others=[2, 0]),
    '265/298': dict(p=265, q=298, Ts=[[-563,-298,-265,265],[-563,-33,265],[-265,265,298,563]], hub=0, others=[1, 2]),
}
if __name__ == '__main__':
    for s in sys.argv[1:]:
        certify(**CONFIGS[s])
