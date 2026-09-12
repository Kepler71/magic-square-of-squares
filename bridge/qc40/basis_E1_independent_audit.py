#!/usr/bin/env python3
# basis_E1_independent_audit.py -- НЕЗАВИСИМАЯ перепроверка выводов basis_E1.sage
# БЕЗ Sage: чистый Python, целые и Fraction, групповой закон выписан руками.
#
# Запуск:  python3 /home/kep/magicKube/bridge/qc40/basis_E1_independent_audit.py
# Выход:   basis_E1_independent_audit.json, basis_E1_independent_audit.log
#
# Что проверяется независимо от Sage/eclib/PARI:
#   [1] P и T лежат на E1 (точная целочисленная подстановка, невязка 0)      -> proved
#   [2] phi1-образ известной точки D (0, 2023) равен РОВНО -2P               -> proved
#       (контроль ловушки индекса 2: образующая -- P, а не образ точки D)
#   [3] E1(Q)_tors = Z/2: gcd(#E(F_p)) по хорошим нечётным p равен 2,
#       и T -- точка порядка 2. Кручение инъективно вкладывается в E(F_p).   -> proved
#   [4] для q = 2,3,5,7 ни P, ни P+T не q-делимы в E1(Q): образ вне q*E(F_p),
#       группа E(F_p) и множество q*E(F_p) построены полным перебором        -> proved
#   [5] hhat(P) = lim h(x(2^k P))/4^k -- третья реализация высоты            -> numerical
#   [6] поиск точек малой высоты в узком окне                                -> observation
#
# ЧЕГО ЭТОТ ФАЙЛ НЕ ДОКАЗЫВАЕТ:
#   - конечности списка простых {2,3,5,7} (граница индекса n <= 8 берётся из
#     Cremona-Siksek в Sage и, независимо, из CPS-границы eclib: см.
#     basis_E1_eclib_bound.log, где eclib сама берёт простые {2,3});
#   - ранга E1(Q) = 1 (это 2-спуск, зависимость от ПО остаётся).
from fractions import Fraction as F
from math import gcd, log
from pathlib import Path
import json

OUT = Path('/home/kep/magicKube/bridge/qc40')
LOG = []


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print(s, flush=True)


# минимальная модель E1: y^2 = x^3 + A2 x^2 + A4 x + A6  (a1 = a3 = 0)
A2, A4, A6 = -1, 1617792974488319, -14262298467262390622975
# сырая модель (new_bridge_verified.json, elliptic_factors[0].raw_ainvs)
R2, R4, R6 = 37608911, 2089269703356959, 7989037884942063852049
a0, a2c, a4c, a6c = 4092529, 47287151, 37608911, 44182609
SHIFT = 12536304                      # X_min = X_raw + SHIFT
Px, Py = 34384993, 286391046720       # образующая свободной части
Tx, Ty = 8443775, 0                   # образующая кручения

out = {'model': {'E1_min_ainvs': [0, A2, 0, A4, A6], 'E1_raw_ainvs': [0, R2, 0, R4, R6],
                 'shift_Xmin_minus_Xraw': SHIFT, 'P': [Px, Py], 'T': [Tx, Ty]}}


def resid(x, y, a2, a4, a6):
    return y * y - (x * x * x + a2 * x * x + a4 * x + a6)


# ---------------------------------------------- [1] точная подстановка
rP = resid(Px, Py, A2, A4, A6)
rT = resid(Tx, Ty, A2, A4, A6)
rPraw = resid(Px - SHIFT, Py, R2, R4, R6)
say("[1] невязка P (min) = %s ; T (min) = %s ; P (raw) = %s  [proved]" % (rP, rT, rPraw))
assert rP == 0 and rT == 0 and rPraw == 0
out['exact_substitution'] = {'P_min': rP, 'T_min': rT, 'P_raw': rPraw, 'status': 'proved'}

# ------------------------------- групповой закон над Q (Fraction), руками
def add(Pp, Q, a2=A2, a4=A4):
    if Pp is None:
        return Q
    if Q is None:
        return Pp
    x1, y1 = Pp
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if Pp == Q:
        if y1 == 0:
            return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - a2 - x1 - x2
    return (x3, lam * (x1 - x3) - y1)


def mul(n, Pp):
    if n < 0:
        Pp = None if Pp is None else (Pp[0], -Pp[1])
        n = -n
    Rr, Q = None, Pp
    while n:
        if n & 1:
            Rr = add(Rr, Q)
        Q = add(Q, Q)
        n >>= 1
    return Rr


PP, TT = (F(Px), F(Py)), (F(Tx), F(Ty))
assert add(TT, TT) is None
P2 = mul(2, PP)
PT = add(PP, TT)
say("[1] 2P = (%s, %s) ; P+T = (%s, %s)  [proved]" % (P2[0], P2[1], PT[0], PT[1]))

# -------------------- [2] КОНТРОЛЬ: образ известной точки D под phi1
# D(Q) ∋ (0, ±2023) при t = 1 -- вырожденная девятка (u0,u4,u8) = (17,17,17).
# phi1 = (a6 x^2, a6 y) -> (0, ±a6*2023) на сырой модели -> X_min = SHIFT.
KX, KY = F(SHIFT), F(a6c * 2023)
assert resid(KX, KY, A2, A4, A6) == 0
coords = None
for k in range(-6, 7):
    for tn, tt in [('O', None), ('T', TT)]:
        Q = mul(k, PP) if tt is None else add(mul(k, PP), tt)
        if Q is not None and (Q[0], Q[1]) == (KX, KY):
            coords = (k, tn)
say("[2] phi1(0, 2023) = (%s, %s) на min-модели; в базисе <P>+<T> это %s  [proved]" % (KX, KY, coords))
assert coords == (-2, 'O'), "образ известной точки НЕ равен -2P -- пересчитать!"
say("[2] КОНТРОЛЬ ЛОВУШКИ: образ известной точки D равен -2P, то есть лежит в 2*E1(Q).")
say("[2] Если бы его взяли за образующую, индекс был бы 2 и решето убивало бы кандидатов "
    "НЕПРАВОМЕРНО. Образующая -- P.  [proved]")
out['known_point_control'] = {'phi1_image_min': [str(KX), str(KY)],
                              'basis_coordinates': {'k': -2, 'torsion': 'O'},
                              'trap': 'образ известной точки = -2P; как образующую использовать НЕЛЬЗЯ',
                              'status': 'proved'}

# ---------------------------------------- дискриминант минимальной модели
b2, b4, b6 = 4 * A2, 2 * A4, 4 * A6
b8 = b2 * A6 - A4 * A4
Delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
DISC_FACT = -(2 ** 20) * (3 ** 2) * (5 ** 2) * (7 ** 2) * (13 ** 2) * (17 ** 8) * (23 ** 6) \
            * (37 ** 2) * (1213 ** 2) * (9397 ** 2)
say("[1] Delta(min) = %s ; совпал с факторизацией лога: %s  [proved]" % (Delta, Delta == DISC_FACT))
assert Delta == DISC_FACT
out['discriminant'] = {'value': str(Delta), 'matches_log_factorization': True, 'status': 'proved'}


# -------------------------------------------- арифметика над F_p (перебор)
def points_mod_p(p):
    a2, a4, a6 = A2 % p, A4 % p, A6 % p
    sq = {}
    for y in range(p):
        sq.setdefault((y * y) % p, []).append(y)
    pts = [None]
    for x in range(p):
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % p
        for y in sq.get(rhs, []):
            pts.append((x, y))
    return pts


def addp(Pp, Q, p):
    a2, a4 = A2 % p, A4 % p
    if Pp is None:
        return Q
    if Q is None:
        return Pp
    x1, y1 = Pp
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if Pp == Q:
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * pow(2 * y1 % p, -1, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
    x3 = (lam * lam - a2 - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def mulp(k, Pp, p):
    if Pp is None:
        return None
    Rr, Q = None, Pp
    while k:
        if k & 1:
            Rr = addp(Rr, Q, p)
        Q = addp(Q, Q, p)
        k >>= 1
    return Rr


# ------------------------- [3] кручение: gcd порядков E(F_p) для хороших p
orders = {}
g = 0
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73]:
    if Delta % p == 0:
        continue
    n = len(points_mod_p(p))
    orders[p] = n
    g = gcd(g, n)
say("[3] #E1(F_p) = %s" % orders)
say("[3] gcd = %s => |E1(Q)_tors| делит %s; T имеет порядок 2 => E1(Q)_tors = Z/2  [proved]" % (g, g))
assert g == 2
say("[3] P != O, P != T (y(P) != 0) => P неторсионна => rank E1(Q) >= 1  [proved]")
out['torsion'] = {'orders_mod_p': orders, 'gcd': g, 'conclusion': 'Z/2',
                  'argument': 'кручение инъективно в E(F_p) при хорошей редукции и p нечётном',
                  'status': 'proved'}

# ---------- [4] независимые сертификаты q-неделимости (полный перебор F_p)
targets = {'P': PP, 'P+T': PT}
certs = {}
for q in [2, 3, 5, 7]:
    certs[str(q)] = {}
    for name, Q in targets.items():
        found = None
        for p in [11, 13, 19, 29, 31, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
            if Delta % p == 0 or p == q:
                continue
            xn, xd = Q[0].numerator, Q[0].denominator
            yn, yd = Q[1].numerator, Q[1].denominator
            if xd % p == 0 or yd % p == 0:
                continue
            xp, yp = xn * pow(xd, -1, p) % p, yn * pow(yd, -1, p) % p
            assert (yp * yp - (xp ** 3 + (A2 % p) * xp * xp + (A4 % p) * xp + A6 % p)) % p == 0
            pts = points_mod_p(p)
            qE = {mulp(q, Rr, p) for Rr in pts}
            if (xp, yp) not in qE:
                found = {'p': p, 'order_E_Fp': len(pts), 'image': [xp, yp], 'size_qE': len(qE),
                         'reason': 'образ %s не лежит в %d*E(F_%d)' % (name, q, p)}
                break
        certs[str(q)][name] = found
        say("[4] q=%d, %s : %s  [proved]" % (q, name, found))
        assert found is not None
say("[4] => ни одно простое q из {2,3,5,7} не делит индекс <P>+<T> в E1(Q)  [proved]")
out['independent_nondivisibility'] = {'certificates': certs,
                                      'note': 'если Q = q*S в E(Q), то образ Q лежит в q*E(F_p) для любого '
                                              'p хорошей редукции; обратное направление не используется',
                                      'status': 'proved'}

# ------------------------------- [5] третья реализация канонической высоты
seq, Q = [], PP
for k in range(0, 9):
    xn, xd = Q[0].numerator, Q[0].denominator
    seq.append(log(max(abs(xn), abs(xd))) / 4 ** k)
    Q = add(Q, Q)
say("[5] h(x(2^k P))/4^k, k=0..8 : %s" % [round(z, 9) for z in seq])
say("[5] предел ~ %.9f ; Sage и PARI дают 1.94633297014005  [numerical]" % seq[-1])
out['canonical_height'] = {'sequence': seq, 'limit_estimate': seq[-1],
                          'sage_pari_value': 1.94633297014005, 'status': 'numerical'}

# ----------------------------------------------- [6] поиск точек малой высоты
small = []
for v in range(1, 60):
    v2 = v * v
    for u in range(-4000, 4001):
        if gcd(u, v) != 1:
            continue
        rhs = u ** 3 + A2 * u * u * v2 + A4 * u * v2 * v2 + A6 * v2 ** 3
        if rhs < 0:
            continue
        w = int(rhs ** 0.5)
        for cand in (w - 1, w, w + 1):
            if cand >= 0 and cand * cand == rhs:
                small.append((u, v, cand))
say("[6] точки с x = u/v^2, |u| <= 4000, v <= 59 : %s  [наблюдение, НЕ доказательство]" % small)
out['small_point_search'] = {'window': '|u|<=4000, v<=59', 'found': [list(s) for s in small],
                             'status': 'observation'}

say("")
say("ИТОГ НЕЗАВИСИМОГО АУДИТА:")
say("  подтверждено без Sage: P, T на кривой; кручение = Z/2; P неторсионна;")
say("  phi1(известная точка D) = -2P; неделимость P и P+T на 2,3,5,7.")
say("  НЕ подтверждено независимо: ранг = 1 (2-спуск) и конечность списка простых")
say("  (граница индекса). Для границы есть две независимые реализации:")
say("  Sage Cremona-Siksek (n <= 8) и eclib CPS (простые {2,3}) -- см. basis_E1_eclib_bound.log.")

out['summary'] = {
    'independently_confirmed': ['P, T на кривой (точная арифметика)',
                                'E1(Q)_tors = Z/2', 'P неторсионна => rank >= 1',
                                'phi1(0,2023) = -2P',
                                'P и P+T не делятся на 2,3,5,7'],
    'NOT_independently_confirmed': ['rank E1(Q) = 1 (зависит от 2-спуска в ПО)',
                                    'граница индекса n <= 8 (Cremona-Siksek в Sage; '
                                    'независимо eclib CPS даёт простые {2,3})'],
    'status': 'proved + proved_software'}

(OUT / 'basis_E1_independent_audit.json').write_text(
    json.dumps(out, indent=2, ensure_ascii=False, default=str) + '\n')
(OUT / 'basis_E1_independent_audit.log').write_text("\n".join(LOG) + "\n")
print("записано:", OUT / 'basis_E1_independent_audit.json', OUT / 'basis_E1_independent_audit.log')
