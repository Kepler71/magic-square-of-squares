# Fable 15.09. Численная проверка: Λ(R) := h(u_+(R)) − h(u_+(−R)) = −2⟨R,S⟩_{E₃} + O(1), S = сумма полюсов u_+.
from sage.all import *
import sys, json
sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable')
from pairing_theory import *


def point_on_E3(f3, zv, y3v):
    w = 1 / (zv + 1 / f3.c0)
    return (f3.L * w, f3.L * f3.c0 * w**2 * y3v)


def sum_of_place(P, place_poly_z, y3_expr_coeffs):
    """сумма точек E₃ над местом: z — корни place_poly_z (над числовым полем), y3 = линейная/квадр. функция z."""
    f3 = P['f3']
    pol = place_poly_z
    if pol.degree() == 1:
        K = QQ; roots = [-pol[0] / pol[1]]
    else:
        K = NumberField(pol, 'a'); roots = [rt for rt, m in pol.change_ring(K).roots()]
        assert len(roots) == pol.degree(), 'поле не расщепляет — берём поле разложения'
    EK = f3.E.change_ring(K)
    tot = EK(0)
    for zv in roots:
        y3v = sum(c * zv**i for i, c in enumerate(y3_expr_coeffs))
        X, Y = point_on_E3(f3, zv, y3v)
        tot += EK(X, Y)
    return tot


def eval_u_conic(P, tau, sign=+1):
    D = P['D']; b = abs(D[0])
    if len(D) == 2:
        zv = 2 * tau / (b * (1 + tau**2)); y3 = (1 - tau**2) / (1 + tau**2)
    else:
        zv = (tau**2 - 1) / D[0]; y3 = tau
    try:
        Av = P['A'].numerator()(zv) / P['A'].denominator()(zv); Bv = P['B'].numerator()(zv) / P['B'].denominator()(zv)
    except ZeroDivisionError:
        return None
    return Av + sign * Bv * y3, zv


def analyze_genus0(Ta, Tb, hmax=60, verbose=True):
    P = setup_pair(Ta, Tb); res = []
    import random; random.seed(3)
    taus = [QQ(p) / q for p in range(-15, 16) for q in range(1, 16) if gcd(p, q) == 1] + \
           [QQ(random.randint(-10**e, 10**e)) / random.randint(1, 10**e) for e in range(2, 13) for _ in range(6)]
    for tau in taus:
        r = eval_u_conic(P, tau, +1); r2 = eval_u_conic(P, tau, -1)
        if r is None or r2 is None or r[0] == 0 or r2[0] == 0: continue
        up, zv = r; um = r2[0]
        res.append(dict(htau=naive_h(tau), hz=naive_h(zv), Lam=naive_h(up) - naive_h(um), hup=naive_h(up)))
    if res and verbose:
        print(f'  род 0: точек {len(res)}; |Λ| до {max(abs(d["Lam"]) for d in res):.2f}; h(u_+) − 2h(z): '
              f'{min(d["hup"]-2*d["hz"] for d in res):.2f}…{max(d["hup"]-2*d["hz"] for d in res):.2f}')
    return dict(Ta=[int(t) for t in Ta], Tb=[int(t) for t in Tb], D=[int(t) for t in P['D']], genus=0, samples=res)


def analyze(Ta, Tb, nmax=12, verbose=True):
    P = setup_pair(Ta, Tb); f3 = P['f3']
    u, Dp, Dz, L = pole_sum(P)
    # места полюсов: Sage даёт place через (многочлен по z, y3 − g(z)); разберём вручную
    S = None
    E3 = f3.E
    tot = E3(0)
    zz = L.base_field().gen(); y3f = L.gen()
    for place, mult in Dp.list():
        if place.is_infinite_place():     # z = ∞ ⇔ X = 0 на модели weier: точки (0,0) (|D|=3) или пара над X=0 (|D|=4)
            # вклад: для |D|=3 — 2-кручение (0,0) с кратностью; для |D|=4 — пара ±, сумма O
            if len(P['D']) == 3: tot += mult * E3(0, 0)
            continue
        k, fr, to = place.residue_field()
        Xk, Yk = to(zz), to(y3f)          # z, y3 в поле вычетов
        if k is QQ or k.degree() == 1:
            zv, y3v = QQ(Xk), QQ(Yk); X, Y = point_on_E3(f3, zv, y3v); pt = E3(X, Y); tot += mult * pt; continue
        Kg = k.galois_closure('b') if k.degree() > 2 else k
        EK = E3.change_ring(Kg)
        embs = k.embeddings(Kg)
        for e in embs:
            zv, y3v = e(Xk), e(Yk); X, Y = point_on_E3(f3, zv, y3v); tot_k = EK(X, Y)
            tot = E3(0) if False else tot
            # суммируем в EK, потом спускаем
            if 'acc' not in dir(): acc = EK(0)
            acc += mult * tot_k
        tot = tot + E3([QQ(c) for c in acc]) if not acc.is_zero() else tot
        del acc
    S = tot
    hS = float(S.height()); tor = S.has_finite_order()
    if verbose: print(f'  S = сумма полюсов = {S}, кручение: {tor}, ĥ(S) = {hS:.4f}')
    # точки R ∈ E₃(ℚ)
    rk = E3.pari_curve().ellrank(); gens = [E3(list(p)) for p in rk[3]]
    if gens:
        sat, idx, reg = E3.saturation(gens); gens = list(sat)
    tors = E3.torsion_points()
    res = []
    if gens:
        G = gens[0]
        for n in range(-nmax, nmax + 1):
            for t in tors:
                R = n * G + t
                up = eval_u(P, R, +1); um = eval_u(P, R, -1)
                if up is None or um is None or up == 0 or um == 0 or R.is_zero(): continue
                Lam = naive_h(up) - naive_h(um)
                pair = float((R + S).height() - R.height() - S.height()) / 2
                res.append(dict(n=n, hR=float(R.height()), Lam=Lam, RS=pair, err=Lam + 2 * pair,
                                hup=naive_h(up), pred=2 * float(R.height()) - pair))
    if res:
        errs = [d['err'] for d in res]
        if verbose:
            print(f'  ранг E₃ = {len(gens)}, точек {len(res)}; Λ+2⟨R,S⟩: min {min(errs):.2f} max {max(errs):.2f}; '
                  f'|Λ| до {max(abs(d["Lam"]) for d in res):.1f}; h(u_+) − (2ĥ(R) − ⟨R,S⟩): '
                  f'{min(d["hup"]-d["pred"] for d in res):.2f}…{max(d["hup"]-d["pred"] for d in res):.2f}')
    return dict(Ta=[int(t) for t in Ta], Tb=[int(t) for t in Tb], D=[int(t) for t in P['D']], S=str(S), S_torsion=tor, hS=hS, rankE3=len(gens), samples=res)


if __name__ == '__main__':
    r, s = int(sys.argv[1]), int(sys.argv[2])
    out = []
    for ci in range(int(sys.argv[3]) if len(sys.argv) > 3 else 0, int(sys.argv[4]) if len(sys.argv) > 4 else 11):
        cl = CLASSES[BIG[ci]]; Ts = [cells_of(S, r, s) for S in cl]
        pairs = [(i, j) for i in range(len(Ts)) for j in range(i + 1, len(Ts)) if len(set(Ts[i]) ^ set(Ts[j])) in (1, 2, 3, 4)]
        print('класс', ci, 'кратность', len(cl), 'пар рода 1:', len(pairs), flush=True)
        for i, j in pairs:
            print(' пара', Ts[i], Ts[j], 'Δ =', sorted(set(Ts[i]) ^ set(Ts[j])), flush=True)
            try:
                out.append(analyze(Ts[i], Ts[j]) if len(set(Ts[i]) ^ set(Ts[j])) >= 3 else analyze_genus0(Ts[i], Ts[j]))
            except Exception as e:
                print('  ошибка:', repr(e)[:200], flush=True)
    json.dump(out, open(f'/home/kep/magicKube/isogeny_mult/fable/pairing_numeric_{r}_{s}.json', 'w'))
