# Приём Фишера (письмо 22.09.2026): в формуле спаривания точка (1:0) ничем не выделена — можно взять
# ЛЮБУЮ точку (x0:z0) из P^1(k) с g(x0,z0) != 0 (эквивалентно: заменить квартику на эквивалентную
# перед вычислением). Простые делители этого значения входят в множество мест по Замечанию 3.3,
# поэтому выбираем точку, у которой норма значения ЛЕГКО факторизуется.
# Claude, 24.09.2026. Контроль: значение спаривания не зависит от выбора точки (проверяется здесь же).

def _norm_size(F, c):
    nm = QQ(c) if _is_QQ(F.k) else QQ(F.k(c).norm())
    return abs(nm.numerator()) * abs(nm.denominator())

def _factor_ease(n, limit=10^6, hard_digits=45, alarm_s=20):
    """ (полностью_ли_разложено, цифр_в_остатке) для |n|: пробные деления, затем — для остатка
        умеренного размера — настоящая факторизация PARI с ограничением по времени. """
    n = ZZ(abs(n))
    if n == 0: return (False, 10^9)
    rest = n
    for p in primes(limit):
        if rest == 1: break
        while rest % p == 0: rest //= p
    if rest == 1 or rest.is_prime(proof=False): return (True, 0)
    d = len(str(rest))
    if d <= hard_digits:
        try:
            alarm(alarm_s); rest.factor(proof=False); cancel_alarm()
            return (True, 0)
        except Exception:
            try: cancel_alarm()
            except Exception: pass
    return (False, d)

def eval_points(k, bound=10):
    """ примитивные пары (x0,z0): целые |.| <= bound, gcd = 1 (для k != Q — ещё и с sqrt(d)) """
    pts = []
    for z0 in range(0, bound + 1):
        for x0 in range(-bound, bound + 1):
            if x0 == 0 and z0 == 0: continue
            if gcd(x0, z0) != 1: continue
            pts.append((k(x0), k(z0)))
    if not _is_QQ(k):
        w = k.gen()
        for z0 in range(0, 3):
            for x0 in range(-2, 3):
                for c in (w, -w, 1 + w, 1 - w):
                    pts.append((k(x0) + c, k(z0)))
    return pts

def choose_a2(F, quartics, bound=10, verbose=True, want=1):
    """ Перебор точек: значение g(x0,z0) по модулю квадратов, ранжирование по лёгкости факторизации.
        quartics — список квартик одного класса (у нас g2 и g3: по 3.2(v) годится любая).
        Возвращает список кандидатов [(a, (x0,z0), индекс_квартики, цифр_остатка, цифр_нормы)]. """
    cands = []
    for gi, g in enumerate(quartics):
        for (x0, z0) in eval_points(F.k, bound):
            v = quartic_eval(g, x0, z0)
            if v == 0: continue
            nm = _norm_size(F, v)
            ok, restdig = _factor_ease(nm.numerator() * nm.denominator())
            cands.append((v, (x0, z0), gi, 0 if ok else restdig, len(str(nm))))
    cands.sort(key=lambda t: (t[3], t[4]))
    # приведение по модулю квадратов — только для лучших кандидатов (дорого над числовым полем)
    out = []
    for (a, pt, gi, rd, nd) in cands[:max(want, 1) * 3]:
        best = (a, rd, nd)
        try:
            a2 = F.small_rep(a)[0]
            nm2 = _norm_size(F, a2)
            ok2, rd2 = _factor_ease(nm2.numerator() * nm2.denominator())
            if (0 if ok2 else rd2, len(str(nm2))) < (rd, nd):
                best = (a2, 0 if ok2 else rd2, len(str(nm2)))
        except Exception:
            pass
        out.append((best[0], pt, gi, best[1], best[2]))
    out.sort(key=lambda t: (t[3], t[4]))
    cands = out
    if verbose and cands:
        a, pt, gi, rd, nd = cands[0]
        print(f"      приём Фишера: точка ({pt[0]}:{pt[1]}) квартики g{gi+2}, |N| цифр {nd}, "
              f"{'разложилось полностью' if rd == 0 else f'остаток {rd} цифр'}")
    return cands[:want] if want else cands

def transform_quartic(F, g, M):
    """ g' = det(M)^{-2} * (g o M), M = [[x0,u],[z0,w]]: те же инварианты I, J (I(goM)=det^4 I),
        тот же класс в Sel^2. g'(1,0) = g(x0,z0)/det^2 — то же по модулю квадратов. """
    k = F.k
    R = PolynomialRing(k if not _is_QQ(k) else QQ, ['X', 'Z'])
    X, Z = R.gens()
    (x0, u), (z0, w) = M
    xx = x0 * X + u * Z
    zz = z0 * X + w * Z
    d = x0 * w - u * z0
    assert d != 0, "вырожденная матрица"
    P = g[0]*xx^4 + g[1]*xx^3*zz + g[2]*xx^2*zz^2 + g[3]*xx*zz^3 + g[4]*zz^4
    P = P / d^2
    co = [P.coefficient({X: 4 - t, Z: t}) for t in range(5)]
    return [k(c) if not _is_QQ(k) else QQ(c) for c in co]


def complete_matrix(F, x0, z0):
    """ дополняет столбец (x0,z0) до обратимой матрицы над k """
    k = F.k
    for (u, w) in [(0, 1), (1, 0), (1, 1), (-1, 1), (0, -1)]:
        if x0 * k(w) - k(u) * z0 != 0:
            return ((x0, k(u)), (z0, k(w)))
    raise AssertionError("не удалось дополнить столбец")


def pair_trick(F, g1, g2, g3, bound=10, verbose=True, extra_norm=16, ncheck=1):
    """ Спаривание с выбором точки вычисления по приёму Фишера.
        ncheck > 1: пересчёт с другими точками — значения обязаны совпасть (контроль). """
    k = F.k
    F.check_quartic(g1); F.check_quartic(g2); F.check_quartic(g3)
    cands = choose_a2(F, [g2, g3], bound=bound, verbose=verbose, want=max(ncheck, 1) * 4)
    assert cands, "нет ненулевых значений квартик"
    vals = []; used = []
    seen = set()
    for (a2raw, pt, gi, rd, nd) in cands:
        keyf = (str(pt), gi)
        if keyf in seen: continue
        seen.add(keyf)
        try:
            # КОРРЕКТНЫЙ приём: заменяем саму квартику на эквивалентную g -> det^-2 (g o M),
            # где первый столбец M — выбранная точка; gamma1 пересчитывается по новой квартике.
            M = complete_matrix(F, pt[0], pt[1])
            if gi == 0:
                g2n, g3n = transform_quartic(F, g2, M), g3
            else:
                g2n, g3n = g2, transform_quartic(F, g3, M)
            F.check_quartic(g2n); F.check_quartic(g3n)
            gsrc = g2n if gi == 0 else g3n
            a2 = gsrc[0]
            assert a2 != 0, "нулевой старший коэффициент после преобразования"
            gam, m = F.gamma1(g1, g2n, g3n)
            den = lcm([QQ(t).denominator() for c in gam for t in ([c] if _is_QQ(k) else list(k(c)))])
            gam = [c * den for c in gam]
            num = gcd([ZZ(t) for c in gam for t in ([QQ(c)] if _is_QQ(k) else list(k(c))) if t != 0])
            if num > 1:
                gam = [c / num for c in gam]
            pls = F.places_for(g1, gam, a2, extra_norm=extra_norm)
        except Exception as ex:
            if verbose: print(f"      точка {pt}: места не построены ({ex})"); continue
        tot = 0
        for pl in pls:
            x, z = F.local_point(g1, gam, pl)
            gv = gam[0] * x^2 + gam[1] * x * z + gam[2] * z^2
            tot += 1 if hilb(k, a2, gv, pl) == -1 else 0
        vals.append(tot % 2); used.append((pt, gi, len(pls), nd))
        if verbose: print(f"      точка ({pt[0]}:{pt[1]}) g{gi+2}: мест {len(pls)}, спаривание {tot % 2}")
        if len(vals) >= ncheck: break
    assert vals, "ни одна точка не дала множества мест"
    assert all(v == vals[0] for v in vals), f"РАСХОЖДЕНИЕ по точкам: {vals} {used}"
    # замечание Codex 25.09: контроль считается выполненным, только если ДЕЙСТВИТЕЛЬНО получено
    # ncheck независимых значений; раньше при единственном успешном значении assert проходил молча.
    assert len(vals) >= ncheck, (f"контроль неполон: запрошено {ncheck} точек, успешных {len(vals)} "
                                 f"({used}) — результат не считается проверенным")
    return vals[0], used

def install_trick(F, bound=10, ncheck=1, verbose=False):
    """ Подменяет F.pair на версию с приёмом Фишера (выбор точки вычисления).
        Возвращает исходный метод — для сверки старым способом. """
    orig = F.pair
    def _pair(g1, g2, g3, reps=1, verbose=verbose, extra_norm=16):
        v, used = pair_trick(F, g1, g2, g3, bound=bound, verbose=verbose,
                             extra_norm=extra_norm, ncheck=ncheck)
        return v, len(used)
    F.pair = _pair
    F._pair_orig = orig
    return orig
