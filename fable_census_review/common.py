# Fable, 14.09.2026. Независимая проверка переписи s<=200. Код написан до чтения
# скриптов Claude (census.py, g1census.py, kolyvagin.py); формат jsonl взят из самих файлов.
from sage.all import (QQ, ZZ, PolynomialRing, EllipticCurve, pari, alarm, cancel_alarm,
                      AlarmInterrupt, gcd, is_square)
import time, json, sys

R = PolynomialRing(QQ, 'z'); zz = R.gen()


def coeffs(r, s):
    """значения c, входящие в клетки 1 +- c z: r, s, s-r, s+r"""
    return [r, s, s - r, s + r]


def nine_cells(r, s, z):
    z = QQ(z)
    out = [QQ(1)]
    for c in coeffs(r, s):
        out += [1 + c * z, 1 - c * z]
    return out


def classify(r, s, z):
    """что происходит с девятью клетками при данном z"""
    z = QQ(z)
    if z == 0:
        return 'z=0 (все клетки 1)'
    cells = nine_cells(r, s, z)
    sq = [QQ(x).is_square() for x in cells]
    if all(sq):
        return 'ALERT: все девять квадраты'
    zero = [x for x in cells if x == 0]
    if zero:
        return 'клетка 0 (противоположная = 2), квадратов %d/9' % sum(sq)
    return 'квадратов %d/9' % sum(sq)


def with_timeout(f, secs, *a, **k):
    try:
        alarm(secs)
        v = f(*a, **k)
        cancel_alarm()
        return v
    except AlarmInterrupt:
        return 'TIMEOUT'
    except Exception as e:
        cancel_alarm()
        return 'ERR: %s' % (str(e)[:200])


def prove_rank0(E, tL=900, tM=300):
    """Два независимых от PARI ellrank способа доказать ранг 0.
    L_ratio: L(E,1)/Omega точно (Sage, модулярные символы через at1 с оценкой ошибки);
             L != 0  =>  ранг 0 (Колывагин).  Требуем ratio >= 1/2, чтобы ненулевость
             следовала из одной оценки ошибки, без гипотезы о константе Манина.
    mwrank:  верхняя граница ранга через 2-спуск (Cremona), только как второе свидетельство."""
    res = {}
    t0 = time.time()
    v = with_timeout(lambda: E.lseries().L_ratio(), tL)
    res['L_ratio'] = str(v)
    res['L_time'] = round(time.time() - t0, 1)
    res['rank0_by_L'] = (not isinstance(v, str)) and v >= QQ(1) / 2
    t0 = time.time()
    def mw():
        C = E.mwrank_curve()
        C.two_descent(verbose=False, second_descent=True)
        return int(C.rank_bound())
    v = with_timeout(mw, tM)
    res['mwrank_bound'] = str(v)
    res['mwrank_time'] = round(time.time() - t0, 1)
    res['rank0_by_mwrank'] = (not isinstance(v, str)) and v == 0
    return res


def torsion_order(E):
    return int(E.torsion_subgroup().order())


def points_on_y2_P(P, H):
    """все аффинные точки y^2 = P(z) наивной высоты < H (PARI hyperellratpoints), плюс
    число точек на бесконечности, вычисленное отдельно"""
    pts = pari.hyperellratpoints(pari(P), H)
    aff = sorted(set(QQ(p[0]) for p in pts))
    d = P.degree()
    lc = P.leading_coefficient()
    if d == 3:
        ninf = 1
    elif d == 4:
        ninf = 2 if QQ(lc).is_square() else 0
    else:
        raise ValueError
    naff = 0
    for x in aff:
        naff += 1 if P(x) == 0 else 2
    return aff, naff, ninf
