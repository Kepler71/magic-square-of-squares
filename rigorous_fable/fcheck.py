# Fable, 25.09.2026. Арифметическая перепроверка closure_<p>_<q>.json (без mwrank): лемма ранга 2, H0, ящик, пустота списков.
import sys, json
sys.path.insert(0, '/home/kep/magicKube/rigorous_fable')
from fmodels import *

def check(fn):
    D = json.load(open(fn)); C = json.load(open(D['certificate']))
    assert D['gens_certain'] and D['full_saturation_index'] == 1
    E = EllipticCurve(QQ, list(map(QQ, D['E']))); assert D['E'] == C['E']
    gens = [E([QQ(c) for c in P]) for P in D['gens']]; Tor = [E([QQ(c) for c in P]) for P in D['torsion']]
    assert len(gens) == 2 and len(Tor) == E.torsion_order()
    I = RealIntervalField(200); a, b, c = [I(QQ(lo), QQ(hi)) for lo, hi in D['gram_intervals']]
    lam = QQ(D['lambda_lower']); assert lam > 0 and ((a*b - c*c)/(a + b)).lower() > I(lam).upper()
    # интервалы Грама пересчитываем заново из точек (без mwrank)
    Bv = ZZ(D['B']); assert Bv == ZZ(C['B']); err = I(Bv).log()/3
    def hint(P, n=6):
        Q = (2**n)*P
        if Q.is_zero(): return I(0)
        return I(max(abs(Q[0].numerator()), Q[0].denominator())).log()/4**n + I(-1, 1)*err/4**n
    h1, h2 = hint(gens[0]), hint(gens[1]); ip = (hint(gens[0] + gens[1]) - h1 - h2)/2
    assert h1.lower() >= a.lower() and h1.upper() <= a.upper() and h2.lower() >= b.lower() and h2.upper() <= b.upper() and ip.lower() >= c.lower() and ip.upper() <= c.upper()
    # константы высот из мёбиусовых матриц сертификата
    Cs = {str(j): ZZ((err + I(ZZ(rec['mobius_norm'])).log()).upper().ceil()) for j, rec in enumerate(C['images'])}
    assert {k: int(v) for k, v in Cs.items()} == D['height_constants']
    hub = D['hub']; j2, j0 = D['others']; cp = {str(p['j']): p['integer_bound'] for p in C['pairs']}
    assert D['pairing_constants'] == {str(j2): cp[str(j2)], str(j0): cp[str(j0)]}
    c2, c0 = QQ(cp[str(j2)]), QQ(cp[str(j0)]); C1, C2, C0 = Cs[str(hub)], Cs[str(j2)], Cs[str(j0)]
    def kval(H):
        am = max(c2, c0)**2/(H - C1)
        if H - max(C0, C2) - am <= 0: return None
        return (c2 + c0)**2/(H - C1) + (C0 + C2 + am)**2/(4*(H - max(C0, C2) - am))
    H0 = ZZ(D['H0']); assert H0 > max(Cs.values()) and kval(H0) is not None and kval(H0) < lam and QQ(D['K_H0']) == kval(H0)
    assert kval(H0 - 1) is None or kval(H0 - 1) >= lam   # H0 минимально (не обязательно для доказательства, контроль)
    bound = ZZ(D['box_bound']); assert bound**2*lam > H0 + C1
    assert D['case_A_combinations'] == len(Tor)*(2*bound + 1)**2
    assert D['sieve'][-1]['remaining'] == len(D['exact_candidates'])
    assert not D['case_A_solutions'] and not D['case_B_solutions']
    S = cells(*map(ZZ, D['slope'].split('/')))
    def qualifies(z):
        return z is not None and z != 0 and all(1 + k*z > 0 and (1 + k*z).is_square() for k in S)
    zs = set()
    for rec in D['exact_candidates']:
        z = None if rec['z'] == 'None' else QQ(rec['z']); assert not qualifies(z); zs.add(rec['z'])
    for zs_ in D['case_B_candidates']:
        z = None if zs_ == 'None' else QQ(zs_); assert not qualifies(z)
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
    for ps in D['case_B_polynomials']: assert Rz(sage_eval(ps, locals={'z': z})) != 0
    assert len(D['case_B_polynomials']) == len(Tor)
    other = [w for w in zs if w not in ('0', 'None')]
    print(fn, 'OK: lambda', lam, 'H0', H0, 'box', bound, 'combinations', D['case_A_combinations'], 'exact candidates', len(D['exact_candidates']),
          '(z=0:', sum(1 for r in D['exact_candidates'] if r['z'] == '0'), ', z=inf:', sum(1 for r in D['exact_candidates'] if r['z'] == 'None'),
          ', other finite:', len(other), ') case B candidates', D['case_B_candidates'])

if __name__ == '__main__':
    for fn in sys.argv[1:]: check(fn)
