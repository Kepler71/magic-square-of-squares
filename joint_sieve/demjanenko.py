# Fable, 15.09.2026. Метод Демьяненко–Манина для кривой наклона через симметрию z ↦ −z.
# S = −S — набор клеток (1+λz квадраты для λ ∈ S); T ⊂ S, T ≠ −T, |T| ∈ {3,4}, E_T ранга 1: E_T(ℚ) = ℤP₀ ⊕ Tors.
# Если z — решение, то и −z; P_T(z) = nP₀+t, P_T(−z) = n'P₀+t'. h(z)=h(−z), |ĥ(Q) − h(x(Q))| ≤ B (Силверман, код Sage),
# |h(x(Q)) − h(z(Q))| ≤ c₁ (мёбиусова связь) ⇒ |n²−n'²|·ĥ(P₀) ≤ 2B+2c₁ =: C.
# Либо |n|,|n'| ≤ M₀=(C/ĥ(P₀)+1)/2 (случай A, перебор), либо n' = ±n ⇒ P_T(−z) = ±P_T(z)+t'' ⇒ z(P+t'') = −z(P) (случай B, алгебра).
from sage.all import *
import sys, json, time, itertools
sys.path.insert(0, '/home/kep/magicKube/joint_sieve')
from jsieve import Factor, slope_cells


def naive_h(x):
    x = QQ(x); return float(log(max(abs(x.numerator()), abs(x.denominator()))))


def check_normalization():
    E = EllipticCurve('37a1'); P = E.gens()[0]; Q = 13 * P
    ratio = naive_h(Q[0]) / 169 / float(P.height())
    assert 0.9 < ratio < 1.1, ratio   # ĥ_Sage ≈ h(x)
    return ratio


def all_square(z, S):
    z = QQ(z)
    return all((1 + lam * z) >= 0 and (1 + lam * z).is_square() for lam in S)


def mobius_z_to_xmin(f):
    """z ↦ X = L c₀/(c₀ z + 1) ↦ x_min = u² X + r (изоморфизм на минимальную модель): матрица с целыми элементами."""
    Emin = f.E.minimal_model(); iso = f.E.isomorphism_to(Emin)
    u, r, s_, t_ = iso.tuple()   # x_min = u^2 x + r ? (Sage: (u,r,s,t) с x = u^2 x' + r) — проверяем численно ниже
    # z = u/v; X = L c0 v/(c0 u + v); x_min = a*X + b (a = u_iso^{-2} или u_iso^2, b = ±r) — определим по образцу
    P = f.gens[0]; xm = iso(P)[0]; X = P[0]
    # ищем a, b с xm = a X + b по двум точкам
    P2 = 2 * P; xm2 = iso(P2)[0]; X2 = P2[0]
    a = (xm - xm2) / (X - X2); b = xm - a * X
    assert iso(3 * P)[0] == a * (3 * P)[0] + b
    # x_min = a L c0 v/(c0 u + v) + b = (a L c0 v + b c0 u + b v)/(c0 u + v): матрица [[b c0, a L c0 + b],[c0, 1]]
    M = matrix(QQ, [[b * f.c0, a * f.L * f.c0 + b], [f.c0, 1]])
    d = lcm([x.denominator() for x in M.list()]); M = M * d
    c1 = float(log(2 * max(abs(x) for x in M.list())))
    return Emin, iso, M, c1


def candidates_caseB(f, S):
    """Рациональные X с z(P+t) = −z(P) для всех t ∈ Tors (включая t=0: z=0). Возвращает список z и флаг «уравнение тривиально»."""
    R = PolynomialRing(QQ, ['X', 'Y']); X, Y = R.gens()
    h = X**3 + f.a2 * X**2 + f.a4 * X + f.a6
    zs = set([QQ(0)]); trivial = []
    L, c0 = QQ(f.L), QQ(f.c0)
    for t in f.tdict:
        if t.is_zero(): continue
        Xt, Yt = QQ(t[0]), QQ(t[1])
        D = (X - Xt)**2
        Np = (h - 2 * Yt * Y + Yt**2) - (f.a2 + X + Xt) * D     # X' D
        cond = c0 * L * (X * D + Np) - 2 * X * Np                 # c0 L (X + X') = 2 X X', умножено на D
        A = cond.coefficient({Y: 0}); B = cond.coefficient({Y: 1})
        assert cond == A + B * Y
        Rx = PolynomialRing(QQ, 'x')
        conv = lambda p: Rx(p.univariate_polynomial().list()) if p != 0 else Rx(0)
        A = conv(A); B = conv(B)
        hx = Rx([f.a6, f.a4, f.a2, 1])
        pol = A if B == 0 else A**2 - B**2 * hx
        if pol == 0:
            trivial.append(str(t)); continue
        for (root, m) in pol.roots():
            if root == 0: zs.add(infinity)
            else: zs.add(-1 / c0 + L / root)
        zs.add(-1 / c0)   # P + t = O или X = Xt: точки кручения (покрыты случаем A), добавляем z(O) для полноты
    return sorted(zs, key=str), trivial


def run(S, T, verbose=True):
    S = sorted(ZZ(x) for x in S); assert sorted(-x for x in S) == S
    f = Factor(T)
    assert set(f.T) <= set(S) and sorted(-x for x in f.T) != f.T
    if not (f.rank_proved and f.rank == 1): return dict(T=list(map(int, f.T)), ok=False, why=f'ранг {f.lo}/{f.hi}, найдено {f.rank}')
    P0 = f.gens[0]; hP0 = float(P0.height())
    Emin, iso, M, c1 = mobius_z_to_xmin(f)
    B = float(Emin.silverman_height_bound())
    SAFETY = 2.0   # запас: граница Силвермана удвоена (константы взяты из кода Sage, статья не открывалась)
    C = 2 * SAFETY * B + 2 * c1; K = C / hP0; M0 = int((K + 1) / 2) + 1
    # случай A
    cands = set()
    for n in range(-M0, M0 + 1):
        Q = n * P0
        for t in f.tdict:
            z = f.z_of(Q + t); cands.add(z)
    # случай B
    zB, trivial = candidates_caseB(f, S)
    cands |= set(zB)
    sols = [z for z in cands if z != infinity and z != 0 and all_square(z, S)]
    res = dict(T=list(map(int, f.T)), ok=True, rank=f.rank, tors=f.tinv, sat_index=int(f.sat_index), hP0=hP0, B=B, c1=c1, C=C, K=K, M0=M0,
               nA=(2 * M0 + 1) * len(f.tdict), nB=len(zB), trivial=trivial, ncand=len(cands),
               nondeg_solutions=[str(z) for z in sols])
    if verbose:
        print(f"  T={f.name}: ĥ(P₀)={hP0:.3f}, B={B:.2f}, c1={c1:.2f}, K={K:.1f}, M0={M0}; кандидатов A={res['nA']}, B={len(zB)}"
              f"{' (ТРИВИАЛЬНО для t=' + str(trivial) + ')' if trivial else ''}; невырожденных решений: {sols}", flush=True)
    return res


def run_slope(r, s, nT=3, verbose=True):
    S = slope_cells(r, s)
    out = []
    for k in (3, 4):
        for T in itertools.combinations(sorted(S), k):
            if sorted(-x for x in T) == list(T): continue
            try:
                f = Factor(T)
            except Exception as e:
                continue
            if not (f.rank_proved and f.rank == 1): continue
            res = run(S, T, verbose=verbose); out.append(res)
            if len([o for o in out if o['ok'] and not o['trivial']]) >= nT: return out
    return out


if __name__ == '__main__':
    print('нормировка ĥ_Sage/h(x):', check_normalization())
    if sys.argv[1] == 'sallows':
        # положительный контроль: клетки ±825, ±578 (из квадрата Саллоуса, p=138600, p−q=97104, gcd 168), z = 168/425²
        S = [825, -825, 578, -578]; z = QQ(168) / 425**2; assert all_square(z, S)
        for T in ([-825, -578, 578], [-825, -578, 825], [-578, 578, 825], [-825, 578, 825]):
            res = run(S, T)
            if res['ok']: print('   известная точка найдена:', str(z) in res['nondeg_solutions'])
    else:
        r, s = int(sys.argv[1]), int(sys.argv[2])
        t0 = time.time(); out = run_slope(r, s, nT=int(sys.argv[3]) if len(sys.argv) > 3 else 3)
        good = [o for o in out if o['ok'] and not o['trivial']]
        closed = bool(good) and all(o['nondeg_solutions'] == [] for o in good)
        print(f'ИТОГ {r}/{s}: множителей ранга 1 использовано {len(good)}; ЗАКРЫТ (нет невырожденных решений): {closed}; время {time.time()-t0:.1f} с', flush=True)
        json.dump(dict(slope=f'{r}/{s}', results=out, closed=closed), open(f'dem_{r}_{s}.json', 'w'))
