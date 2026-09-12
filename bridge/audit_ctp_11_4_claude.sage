# НЕЗАВИСИМЫЙ ПЕРЕСЧЁТ спаривания Касселса-Тейта для G1 (m,n)=(11,4).
# Ничего не читает из существующих сертификатов: кривая, класс delta и ВСЕ ТРИ квартики
# строятся заново из (m,n) через descent/ctp_quartic.sage (quartic_from_delta).
# Claude, аудит 2026-09-12.

import json, time, functools
print = functools.partial(print, flush=True)

load('/home/kep/magicKube/descent/ctp_quartic.sage')

T0 = time.time()
def tick(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}")

m, n = 11, 4
s = QQ(m^2 + n^2) / 2
b = s * m^2 * n^2
tick(f"m,n = {m},{n}   s = {s}   b = {b}")

# E : V^2 = (X+b)(X + s m^4)(X + s n^4)
R.<X> = QQ[]
f = (X + b) * (X + s*m^4) * (X + s*n^4)
E0 = EllipticCurve([0, f[2], 0, f[1], f[0]])
roots0 = [-b, -s*m^4, -s*n^4]
for e in roots0:
    assert f(e) == 0
tick(f"E0 = {E0}")

M = E0.minimal_model()
iso = E0.isomorphism_to(M)
tick(f"M  = {M}")
# перенос корней в модель M (порядок ОБЯЗАН сохраниться)
u, r_, s_, t_ = iso.u, iso.r, iso.s, iso.t
print("   iso (u,r,s,t) =", (u, r_, s_, t_))
roots = [(e - r_) / u^2 for e in roots0]
print("   корни M =", roots)
fM = M.division_polynomial(2) / 4
for e in roots:
    assert fM(e) == 0, "корень не перенесён"
scale = u^2
print("   x_E0 = u^2 x_M + r,  u^2 =", scale, " квадрат?", scale.is_square())

I, J = IJ_of_curve(M)
tick(f"I = {I}\n           J = {J}")
F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, e) for e in roots]
print("   phi =", phis)
order = _comp_order(F, phis)
print("   порядок компонент L ->", order)

def mk(vals):
    """vals = (v1,v2,v3) в порядке корней (e1,e2,e3) -> элемент L"""
    return F.E.from_comps([QQ(vals[i]) for i in order])

def cls(z):
    return [QQ(F.E.comp(z, i)).squarefree_part() for i in range(3)]

def cls_in_root_order(z):
    c = [QQ(F.E.comp(z, i)).squarefree_part() for i in range(3)]
    out = [None]*3
    for j, i in enumerate(order):
        out[j] = c[i]
    return out

# --- целевой класс ---
# X - e1 = (mn)^2 u4^2,  X - e2 = s m^2 u0^2,  X - e3 = s n^2 u8^2  => delta = (1, s, s)
delta1 = mk([1, s, s])
tick(f"delta целевой (в порядке корней) = {cls_in_root_order(delta1)}   (ожидается [1, 274, 274])")
assert cls_in_root_order(delta1) == [1, 274, 274]

# --- партнёры по спариванию: два РАЗНЫХ элемента Sel^2 ---
partners = [(2055, 137, 15), (30, 137, 4110)]

results = {}
for pi, pv in enumerate(partners):
    tick(f"=== партнёр {pi}: delta2 = {pv} ===")
    delta2 = mk(pv)
    delta3 = delta1 * delta2
    print("   delta3 класс =", cls_in_root_order(delta3))
    gs = []
    for tag, dl in [('g1', delta1), ('g2', delta2), ('g3', delta3)]:
        t1 = time.time()
        g = F.quartic_from_delta(dl, tag=f"{tag}_{pi}")
        tick(f"   {tag} построена за {time.time()-t1:.1f}s: {g}")
        F.check_quartic(g)
        zc = cls_in_root_order(F.z_inv(g))
        print(f"      z({tag}) класс =", zc, " == delta?", zc == cls_in_root_order(dl))
        assert zc == cls_in_root_order(dl)
        gs.append(g)
    g1, g2, g3 = gs
    t1 = time.time()
    val, npl = F.pair(g1, g2, g3, reps=3)
    tick(f"   <g1,g2>_CT = {val}  ({'-1' if val else '+1'}), мест {npl}, {time.time()-t1:.1f}s")
    t1 = time.time()
    valr, nplr = F.pair(g2, g1, g3, reps=2)
    tick(f"   ОБРАТНО <g2,g1>_CT = {valr} ({'-1' if valr else '+1'}), мест {nplr}, {time.time()-t1:.1f}s")
    t1 = time.time()
    val13, npl13 = F.pair(g1, g3, g2, reps=1)
    tick(f"   <g1,g3>_CT = {val13} (ожидается = <g1,g2> = {val})")
    results[str(pv)] = dict(pair=int(val), places=int(npl), reverse=int(valr),
                            places_rev=int(nplr), pair_g1_g3=int(val13),
                            g1=[str(c) for c in g1], g2=[str(c) for c in g2],
                            g3=[str(c) for c in g3])

# --- контроль: элементы образа E(Q)/2E(Q) обязаны спариваться тривиально ---
tick("=== контроль: точки E(Q) ===")
ctrl = []
T2 = [P for P in M.torsion_points() if P.order() == 2]
pts = list(T2)
try:
    gens = M.gens(proof=False)
    pts += gens
    print("   генераторы MW (не доказано):", gens)
except Exception as ex:
    print("   генераторы не получены:", ex)
pv = partners[0]
delta2 = mk(pv); delta3 = delta1 * delta2
g2c = F.quartic_from_delta(delta2, tag='ctrl_g2')
for P in pts:
    xP = P[0]
    dl = mk([xP - e for e in roots])
    c = cls_in_root_order(dl)
    if c == [1, 1, 1]:
        print(f"   x={xP}: тривиальный класс, пропуск")
        continue
    dl3 = dl * delta2
    try:
        gA = F.quartic_from_delta(dl, tag=f"ctrl_{xP}")
        gB = F.quartic_from_delta(dl3, tag=f"ctrl3_{xP}")
        v, _ = F.pair(gA, g2c, gB, reps=1)
        print(f"   x={xP}  класс={c}  <.,g2> = {v}  (ОБЯЗАН быть 0)")
        ctrl.append(dict(x=str(xP), cls=[str(t) for t in c], pair=int(v)))
    except Exception as ex:
        print(f"   x={xP}: не удалось: {ex}")

out = dict(pair="(11,4)", I=str(I), J=str(J), M=str(M.ainvs()),
           roots=[str(r) for r in roots], phis=[str(p) for p in phis],
           delta_target=[1, 274, 274], results=results, controls=ctrl)
open('/home/kep/magicKube/bridge/audit_ctp_11_4_claude.json', 'w').write(json.dumps(out, indent=1, ensure_ascii=False))
tick("ГОТОВО")
