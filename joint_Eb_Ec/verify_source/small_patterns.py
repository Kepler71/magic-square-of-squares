# Независимая точная проверка случая ранга 1 для малых индексов (без теорем о примитивных делителях).
# На E_N: y^2 = x^3 - N^2 x положим x(G) = N t. Тогда x(kG) = N f_k(t), где f_k -- x-отображение
# умножения на k на E_1: y^2 = x^3 - x (изоморфизм над Q(sqrt N) коммутирует с [k], x-отображение над Q).
# Условие прогрессии x(2 k1 G) + x(2 k3 G) = 2 x(2 k2 G)  <=>  F(t) := f_{2k1}+f_{2k3}-2 f_{2k2} = 0,
# многочлен от t НЕ зависит от N. Любой рациональный корень t не из {0, 1, -1} даёт кривую E_N
# (N = бесквадратная часть t^3 - t) и неторсионную G -- т.е. контрпример. Проверяем все тройки
# различных k1,k2,k3 >= 1 (k2 -- середина прогрессии) с max = M <= MMAX.
import os, sys, time, json, itertools
from sage.all import EllipticCurve, QQ, PolynomialRing
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
E = EllipticCurve(QQ, [-1, 0])
R = PolynomialRing(QQ, 't'); t = R.gen()
F = R.fraction_field()
cache = {}
def f(k):
    if k not in cache:
        phi = E.multiplication_by_m(k, x_only=True)   # рациональная функция от x
        num, den = phi.numerator(), phi.denominator()
        cache[k] = F(R(num.change_ring(QQ).subs(x=t)) if False else R(str(num).replace('x','t'))) / F(R(str(den).replace('x','t')))
    return cache[k]
# КОНТРОЛЬ 1: f_k совпадает с прямой арифметикой точек на E_5, G=(-4,6), t=-4/5
E5 = EllipticCurve(QQ, [-25, 0]); G5 = E5(-4, 6); t5 = QQ(-4)/5
for k in range(1, 13):
    assert (k*G5)[0] == 5 * f(k)(t5), k
print("контроль 1 (f_k против арифметики на E_5, k<=12): ok")
# КОНТРОЛЬ 2 (положительный): многочлен с заведомым рациональным корнем t0 находится
t0r = QQ(7)/3
Pc = (f(10) + f(2) - (f(10)(t0r) + f(2)(t0r))).numerator()
assert t0r in [r for r, _ in Pc.roots(QQ)]
print("контроль 2 (поиск корня t0=7/3 степени", Pc.degree(), "): ok")
# КОНТРОЛЬ 3: на E_1 x(2G) = (x^2+1)^2/(4(x^3-x))
assert f(2) == (t**2+1)**2/(4*(t**3-t))
print("контроль 3 (формула удвоения): ok")
results = []
t0 = time.time()
for M in range(3, MMAX + 1):
    for pair in itertools.combinations(range(1, M), 2):
        ks = sorted(pair + (M,))
        for mid in ks:
            ends = [k for k in ks if k != mid]
            expr = f(2*ends[0]) + f(2*ends[1]) - 2*f(2*mid)
            P = expr.numerator()
            roots = [r for r, _ in P.roots(QQ)]
            bad = [r for r in roots if r not in (0, 1, -1)]
            # дополнительно: исключить корни, где знаменатели обнуляются (там выражение не определено)
            real_bad = []
            for r in bad:
                ok = all(f(2*k).denominator()(r) != 0 for k in ks)
                if ok: real_bad.append(str(r))
            results.append({"M": M, "ks": ks, "mid": mid, "deg": int(P.degree()),
                            "rational_roots": [str(r) for r in roots], "counterexample_t": real_bad})
            print(M, ks, "mid", mid, "deg", P.degree(), "rat.roots", [str(r) for r in roots],
                  "КОНТРПРИМЕР" if real_bad else "ok", "%.1fs" % (time.time() - t0), flush=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "small_patterns_M%d.json" % MMAX)
json.dump(results, open(out, "w"), indent=1, ensure_ascii=False)
print("всего шаблонов:", len(results), "контрпримеров:", sum(1 for r in results if r["counterexample_t"]))
