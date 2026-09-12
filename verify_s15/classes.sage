# П.1: полнота списка descent-классов Q: Y^2 = (109X^2-229)(X^4-1), X = p/q, gcd(p,q)=1.
# f1 = p^2-q^2 = d1 r^2, f2 = p^2+q^2 = d2 s^2, f3 = 109p^2-229q^2 = d3 w^2, d1 d2 d3 = квадрат.
# Кандидаты (вывод по НОД форм): d1 in ±{1,2,3,5,6,10,15,30}, d2 in {1,2,13,26}, d3 = sqf(d1 d2).
# Для каждой тройки: известная рациональная точка ИЛИ доказательство отсутствия точек над R или Q_l.
import functools, itertools
print = functools.partial(print, flush=True)
Rx.<x> = QQ[]
F = [x^2 - 1, x^2 + 1, 109*x^2 - 229]          # f_i(p/q, 1) = f_i(p,q)/q^2 (q=1 карта)
Fh = [lambda p, q: p^2 - q^2, lambda p, q: p^2 + q^2, lambda p, q: 109*p^2 - 229*q^2]

def is_sq_Ql(a, l):
    # a in Q (или Z), a != 0: квадрат ли в Q_l
    a = QQ(a)
    v = a.valuation(l)
    if v % 2: return False
    u = a / l^v
    if l == 2:
        un = (u.numerator() * u.denominator()) % 8     # u = n/d, n,d нечётны; n/d ≡ n*d (mod 8)
        return un == 1
    return kronecker(ZZ(u.numerator() * u.denominator()), l) == 1

def class_ok_at(vals, ds, l):
    # vals: значения форм (0 допускается = квадрат 0)
    return all(v == 0 or is_sq_Ql(v / d, l) for v, d in zip(vals, ds))

def locally_solvable(ds, l, K=None):
    # Строгий тест: есть ли (p:q) in P^1(Q_l) с f_i(p,q)/d_i квадратом (или 0) для всех i.
    # Шар p ≡ a (mod l^k): форма f_i "определена" на шаре, если v(f_i(a)) + e <= k (тогда квадратичный класс
    # f_i постоянен на шаре, т.к. f_i(x) - f_i(a) ≡ 0 mod l^k). Если какая-то определённая форма имеет неверный
    # класс — в шаре точек нет (строго). Если все определены и верны — точка есть. Иначе дробим до глубины K.
    e = 3 if l == 2 else 1
    if K is None: K = 16 if l == 2 else 10
    unknown = []
    def rec(chart, a, k):
        p, q = (a, 1) if chart == 0 else (1, a)
        vals = [fh(p, q) for fh in Fh]
        und = []
        for i, (v, d) in enumerate(zip(vals, ds)):
            if v != 0 and v.valuation(l) + e <= k:
                if not is_sq_Ql(v / d, l):
                    return False                      # определённая форма с неверным классом: точек нет
            else:
                und.append(i)
        if not und:
            return True
        if k >= K:
            unknown.append((chart, a, k, und)); return False
        return any(rec(chart, a + j*l^k, k + 1) for j in range(l))
    if rec(0, 0, 0) or rec(1, 0, 1):
        return True, "point found"
    if not unknown:
        return False, "no points (rigorous)"
    return None, f"UNKNOWN balls: {unknown[:3]}"

def real_ok(ds):
    d1, d2, d3 = ds
    if d2 < 0: return False
    # |X| < 1: f1<0, f3<0 ; |X| > sqrt(229/109): f1>0, f3>0 ; между: f3<0<f1
    # d3 = sqf(d1 d2) имеет знак d1, поэтому интервал 1 < X^2 < 229/109 (f1>0>f3) исключается автоматически
    return d1*d3 > 0

def sqf(n):
    n = ZZ(n); s = sign(n)
    return s * prod(pr^(e % 2) for pr, e in n.abs().factor())

known = {(6,26,39): 5, (-6,26,-39): QQ(1)/5, (5,13,65): QQ(3)/2, (-5,13,-65): QQ(2)/3,
         (15,2,30): QQ(17)/7, (-15,2,-30): QQ(26399)/54721}
for cls, X0 in known.items():
    p0, q0 = X0.numerator(), X0.denominator()
    got = tuple(sqf(fh(p0, q0)) for fh in Fh)
    assert got == cls, (cls, got)
print("known points lie in their classes: OK")

# контроль теста: на классах с известной точкой он обязан давать «разрешимо» везде
primes_to_check = [2, 3, 5, 13, 109, 229] + [l for l in primes(7, 100) if l not in (13,)]
bad_ctrl = [(cls, l) for cls in known for l in primes_to_check if locally_solvable(cls, l)[0] is not True]
print("control (known classes solvable at all checked primes):", "OK" if not bad_ctrl else f"FAIL {bad_ctrl}")
D1 = [s_*d for s_ in (1, -1) for d in (1, 2, 3, 5, 6, 10, 15, 30)]
D2 = [1, 2, 13, 26]
primes_to_check = [2, 3, 5, 13, 109, 229] + [l for l in primes(7, 100) if l not in (13,)]
els = []
for d1, d2 in itertools.product(D1, D2):
    d3 = sqf(d1*d2)
    ds = (d1, d2, d3)
    if ds in known:
        els.append(ds); print(f"{str(ds):16s} known rational point X = {known[ds]}"); continue
    if not real_ok(ds):
        print(f"{str(ds):16s} no real points"); continue
    obstr = None
    for l in primes_to_check:
        ok, why = locally_solvable(ds, l)
        if ok is False:
            obstr = (l, why); break
        if ok is None:
            print(f"   {ds} at l={l}: {why}")
    if obstr:
        print(f"{str(ds):16s} no Q_{obstr[0]}-points ({obstr[1]})")
    else:
        els.append(ds); print(f"{str(ds):16s} *** locally solvable at R and all checked primes, no known point ***")
print("\neverywhere-locally-solvable (among checked places):", els)
print("equals Codex's six classes:", sorted(els) == sorted(known.keys()))
