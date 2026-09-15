from demjanenko import *
import random
random.seed(1)
S = slope_cells(204, 247)
facs = []
for T in itertools.combinations(sorted(S), 3):
    if sorted(-x for x in T) == list(T): continue
    f = Factor(T)
    if f.rank_proved and f.rank == 1: facs.append(f)
    if len(facs) >= 6: break
print('множителей ранга 1 для проверки:', [f.name for f in facs])
# (а) случай B над F_p: корни pol vs прямой перебор
R = PolynomialRing(QQ, ['X', 'Y']); X, Y = R.gens()
for f in facs[:3]:
    L, c0 = QQ(f.L), QQ(f.c0)
    for p in [1009, 2003]:
        if p in f.bad: continue
        Fp = GF(p); El = f.E.change_ring(Fp); pts = El.points()
        def zp(P):
            if P.is_zero(): return -1/Fp(c0)
            if P[0] == 0: return 'inf'
            return -1/Fp(c0) + Fp(L)/P[0]
        for t in f.tdict:
            if t.is_zero(): continue
            tp = El(t)
            direct = set()
            for P in pts:
                if P.is_zero() or P == -tp or P[0] == Fp(t[0]): continue   # исключаем X = X_t
                a, b = zp(P), zp(P + tp)
                if a != 'inf' and b != 'inf' and a + b == 0: direct.add(P[0])
            # многочлен из candidates_caseB (повтор алгебры)
            Xt, Yt = QQ(t[0]), QQ(t[1]); h = X**3 + f.a2*X**2 + f.a4*X + f.a6
            D = (X - Xt)**2; Np = (h - 2*Yt*Y + Yt**2) - (f.a2 + X + Xt)*D
            cond = c0*L*(X*D + Np) - 2*X*Np
            A = cond.coefficient({Y: 0}); B = cond.coefficient({Y: 1})
            Rx = PolynomialRing(QQ, 'x'); conv = lambda q: Rx(q.univariate_polynomial().list()) if q != 0 else Rx(0)
            A = conv(A); B = conv(B); hx = Rx([f.a6, f.a4, f.a2, 1]); pol = A if B == 0 else A**2 - B**2*hx
            polp = pol.change_ring(Fp); roots = set(r for r, m in polp.roots())
            ok = direct <= roots
            print(f'  (а) T={f.name}, p={p}, t порядка {t.order()}: прямых решений {len(direct)}, корней pol mod p {len(roots)}, прямые ⊂ корни: {ok}', flush=True)
            assert ok
# (б) высоты: |ĥ(Q) − h(z(Q))| ≤ B + c1 на кратных
for f in facs[:4]:
    Emin, iso, M, c1 = mobius_z_to_xmin(f); B = float(Emin.silverman_height_bound())
    worst = 0
    for n in range(-7, 8):
        for t in f.tdict:
            Q = n*f.gens[0] + t
            if Q.is_zero(): continue
            z = f.z_of(Q)
            if z == infinity: continue
            d = abs(float(Q.height()) - naive_h(z)); worst = max(worst, d)
    print(f'  (б) T={f.name}: max|ĥ(Q) − h(z(Q))| = {worst:.2f} при границе B + c1 = {B + c1:.2f}: {"OK" if worst <= B + c1 else "НАРУШЕНО"}', flush=True)
    assert worst <= B + c1
# (в) сквозной контроль: z с обеими P_T(±z) рациональными среди nP0 + t, |n| ≤ 4
for f in facs:
    zs = {}
    for n in range(-4, 5):
        for t in f.tdict:
            Q = n*f.gens[0] + t; z = f.z_of(Q)
            if z != infinity: zs.setdefault(z, []).append((n, str(t)[:20]))
    sym = [z for z in zs if z != 0 and -z in zs]
    if sym:
        res = run(S, f.T, verbose=False)
        # метод обязан выдать эти z среди кандидатов (проверяем полный список кандидатов, а не только «решения»)
        cands = set()
        for n in range(-res['M0'], res['M0'] + 1):
            for t in f.tdict: cands.add(f.z_of(n*f.gens[0] + t))
        zB, triv = candidates_caseB(f, S); cands |= set(zB)
        found = all(z in cands for z in sym)
        print(f'  (в) T={f.name}: симметричные z (обе P_T(±z) рациональны): {[(str(z), zs[z], zs[-z]) for z in sym][:3]} — все среди кандидатов метода: {found}', flush=True)
        assert found
print('проверки (а)(б)(в) пройдены')
