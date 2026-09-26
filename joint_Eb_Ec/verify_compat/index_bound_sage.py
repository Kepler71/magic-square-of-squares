# Граница индекса [E(Q) : <gens>+tors] через нижнюю оценку канонической высоты
# (Sage HeightFunction.min — собственный код Sage по Cremona–Prickett–Siksek, не eclib/mwrank)
# и константу Эрмита: R_E >= (lam/gamma_r)^r, gamma_1 = 1, gamma_2^2 = 4/3.
from sage.all import *
import json
data = json.load(open('/home/kep/magicKube/joint_Eb_Ec/codex_unified/unified_theorems_20260925/rank_compatibility.json'))
out = []
for cur in data['curves']:
    n = cur['n']
    E = EllipticCurve(QQ, [0, 0, 0, -n*n, 0])
    P = [E([QQ(c) for c in pt]) for pt in cur['points']]
    r = len(P)
    # сравнение нормировок высоты Sage и PARI
    hs = [p.height(precision=100) for p in P]
    hp = [RealField(100)(pari(E).ellheight(list(p.xy()))) for p in P]
    M = matrix(RealField(100), r, r, lambda i, j: (P[i]+P[j]).height(precision=100) - P[i].height(precision=100) - P[j].height(precision=100))
    M = M/2 + diagonal_matrix([hs[i] for i in range(r)])*0  # off-diagonal pairing
    for i in range(r): M[i, i] = hs[i]
    Reg = M.det()
    # ВАЖНО: HeightFunction.min на неминимальной модели y^2=x^3-n^2x (n=3400, 3366) выдаёт
    # «нижнюю границу», превышающую высоты реальных точек, т.е. неверна там; берём минимальную модель.
    lam_nonmin = E.height_function().min(0.0001, 20)
    Em = E.minimal_model(); iso = E.isomorphism_to(Em)
    Pm = [iso(p) for p in P]
    assert all(abs(Pm[i].height(precision=100) - hs[i]) < 1e-20 for i in range(r))
    lam = Em.height_function().min(0.0001, 20)
    # контроль: lam не больше высот небольших комбинаций образующих
    import itertools
    small = min(sum((c*q for c, q in zip(cs, Pm)), Em(0)).height() for cs in itertools.product(range(-2, 3), repeat=r) if any(cs))
    assert lam <= small, (lam, small)
    if r == 1:
        m = sqrt(Reg/lam)
    else:
        m = 2*sqrt(Reg/3)/lam
    rec = dict(n=n, rank=r, minimal_model=list(map(int,Em.a_invariants())), lambda_on_given_model=float(lam_nonmin), given_model_minimal=bool(E.is_minimal()), min_height_small_combos=float(small), sage_heights=[float(h) for h in hs], pari_heights=[float(h) for h in hp],
               regulator=float(Reg), lambda_lower=float(lam), index_bound=float(m))
    print(rec)
    out.append(rec)
json.dump(out, open('/home/kep/magicKube/joint_Eb_Ec/verify_compat/index_bound.json', 'w'), indent=1)
