# Claude, 25.09: независимый аудит строгих закрытий Fable (73/362, 265/298). Код Fable не импортируется;
# из его JSON берутся только входные данные (E, наборы T, матрицы, константы), всё остальное пересчитывается.
from sage.all import *
import json, sys, time
sl = sys.argv[1]; tag = sl.replace('/', '_')
cert = json.load(open(f'/home/kep/magicKube/rigorous_fable/certificate_{tag}.json'))
clo = json.load(open(f'/home/kep/magicKube/rigorous_fable/closure_{tag}.json'))
r, s = map(int, sl.split('/'))
S = [s, -s, r, -r, s - r, r - s, s + r, -(s + r)]
a1, a2, a3, a4, a6 = [QQ(t) for t in cert['E']]
assert a1 == 0 and a3 == 0
E = EllipticCurve([0, a2, 0, a4, a6])
f = lambda X: X**3 + a2*X**2 + a4*X + a6
fp = lambda X: 3*X**2 + 2*a2*X + a4
Ts = [[int(c) for c in T] for T in cert['Ts']]
Ms = [[[QQ(v) for v in row] for row in im['M']] for im in cert['images']]
hub = int(cert['hub']); o1, o2 = [int(t) for t in cert['others']]
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Kz = Rz.fraction_field()
mob = lambda m, t: (m[0][0]*t + m[0][1]) / (m[1][0]*t + m[1][1])
rep = {}
# (а) отображения
ok_maps = []
for i in range(3):
    ratio = Kz(f(mob(Ms[i], Kz(z)))) / Kz(prod(1 + c*z for c in Ts[i]))
    nu, de = ratio.numerator(), ratio.denominator()
    lc = QQ(nu.leading_coefficient()) / QQ(de.leading_coefficient())
    ok_maps.append(bool(Rz(nu/nu.leading_coefficient()).is_square() and Rz(de/de.leading_coefficient()).is_square() and lc.is_square()))
rep['maps'] = ok_maps; print('(а) отображения:', ok_maps, flush=True)
# (б) ранг другим путём (PARI), насыщение, кручение
rk = pari(E).ellrank(); lo, hi = int(rk[0]), int(rk[1])
G = [E([QQ(t) for t in g[:2]]) for g in clo['gens']]
sat, idx, reg = E.saturation(G, max_prime=0)
tors = E.torsion_points()
rep.update(rank=(lo, hi), sat_index=int(idx), ntors=len(tors))
print(f'(б) PARI ellrank [{lo},{hi}], насыщение {idx}, кручение {len(tors)}', flush=True)
# (в) матрица Грама численно против lambda_lower
Gm = matrix(RR, 2, 2, [[G[0].height(), (( (G[0]+G[1]).height() - (G[0]-G[1]).height()) / 4)],
                       [(( (G[0]+G[1]).height() - (G[0]-G[1]).height()) / 4), G[1].height()]])
lam = QQ(clo['lambda_lower']); lmin = min(Gm.eigenvalues())
rep['lambda_min_num'] = float(lmin); rep['lambda_ok'] = bool(lmin > lam)
print(f'(в) λ_min численно {float(lmin):.4f} > λ Fable {lam} = {float(lam):.4f}: {lmin > lam}', flush=True)
# (г) лемма ранга 2 — своя общая формула
Ch = QQ(clo['height_constants'][str(hub)]); Ca = QQ(clo['height_constants'][str(o1)]); Cb = QQ(clo['height_constants'][str(o2)])
P1 = QQ(clo['pairing_constants'][str(o1)]); P2 = QQ(clo['pairing_constants'][str(o2)])
def K(H):
    H = QQ(H); A = max(P1, P2)**2 / (H - Ch)
    return (P1 + P2)**2 / (H - Ch) + (Ca + Cb + A)**2 / (4 * (H - max(Ca, Cb) - A))
H0 = QQ(clo['H0'])
rep['K_H0'] = str(K(H0)); rep['K_ok'] = bool(K(H0) < lam)
print(f'(г) K(H0={H0}) = {float(K(H0)):.6f} < λ = {float(lam):.6f}: {K(H0) < lam}', flush=True)
R2 = (H0 + Ch) / lam                      # ĥ(Q_hub) < H0 + C_hub, ĥ >= λ (n²+m²)
BOX = int(floor(sqrt(R2))) + 1
print(f'   круг n²+m² < {float(R2):.1f}, радиус {BOX} (у Fable ящик {clo["box_bound"]})', flush=True)
# (д) своё решето по хабу на СВОИХ простых
Mh = Ms[hub]; A_, B_, C_, D_ = Mh[0][0], Mh[0][1], Mh[1][0], Mh[1][1]
cands = [(n, m, ti) for n in range(-BOX, BOX + 1) for m in range(-BOX, BOX + 1) if n*n + m*m < R2 for ti in range(len(tors))]
print('   кандидатов:', len(cands), flush=True)
t0 = time.time()
used = []
for p in primes(53, 400):
    if E.discriminant().valuation(p) != 0 or (A_*D_ - B_*C_).valuation(p) != 0: continue
    if any(QQ(v).denominator() % p == 0 for row in Mh for v in row): continue
    Ep = E.change_ring(GF(p)); F = GF(p)
    g1, g2 = Ep(G[0]), Ep(G[1]); tp = [Ep(t) for t in tors]
    m1 = {n: n * g1 for n in range(-BOX, BOX + 1)}; m2 = {m: m * g2 for m in range(-BOX, BOX + 1)}
    keep = []
    for (n, m, ti) in cands:
        P = m1[n] + m2[m] + tp[ti]
        if P.is_zero(): keep.append((n, m, ti)); continue
        den = F(C_) * P[0] - F(A_)
        if den == 0: keep.append((n, m, ti)); continue
        zz = (F(B_) - F(D_) * P[0]) / den
        if all((1 + F(c) * zz).is_square() for c in S): keep.append((n, m, ti))
    cands = keep; used.append((p, len(cands)))
    print(f'   p={p}: осталось {len(cands)}', flush=True)
    if len(cands) <= 6: break
sols = []
for (n, m, ti) in cands:
    P = n * G[0] + m * G[1] + tors[ti]
    if P.is_zero(): continue
    den = C_ * P[0] - A_
    if den == 0: continue
    zz = (B_ - D_ * P[0]) / den
    if zz != 0 and all((1 + c*zz) > 0 and QQ(1 + c*zz).is_square() for c in S): sols.append(str(zz))
rep['sieve'] = used; rep['caseA_solutions'] = sols
print(f'(д) случай A: выжившие {cands} -> невырожденных решений {sols} ({time.time()-t0:.0f}s)', flush=True)
# (е) случай B: Q_o1 − ε Q_o2 = T  =>  x(Q_o2) = x(Q_o1 − T)
Xa = mob(Ms[o1], Kz(z)); Xb = mob(Ms[o2], Kz(z))
zs = set()
for t in tors:
    eq = (Xb - Xa) if t.is_zero() else (Xb - (t[0] + fp(t[0]) / (Xa - t[0])))
    nu = Kz(eq).numerator()
    if nu == 0: print('   ТОЖДЕСТВЕННЫЙ НОЛЬ при T =', t); continue
    zs |= set(rt for rt, _ in nu.roots(QQ))
badB = [str(zz) for zz in zs if zz != 0 and all((1 + c*zz) > 0 and QQ(1 + c*zz).is_square() for c in S)]
rep['caseB_roots'] = sorted(str(t) for t in zs); rep['caseB_solutions'] = badB
print(f'(е) случай B: корни z {sorted(zs)} -> невырожденных решений {badB}', flush=True)
json.dump(rep, open(f'audit_{tag}.json', 'w'), ensure_ascii=False, indent=1, default=str)
