# Fable, 14.09.2026. Проверка леммы §2 THEOREM_INFINITY_CONSTRAINTS (Codex):
# род (Z/2)^r-накрытия P^1, заданного r линейными множителями с различными нулями.
# (а) независимая формула Кани–Розена: g(C) = sum по непустым T подмножествам {1..r} рода y^2 = prod_{i in T} f_i;
# (б) прямой Sage-род проективных моделей при r = 2, 3, 4 (в том числе кривая Саллоуса a=4420, b=-3360
#     и слой A при p0 = 1/3); (в) тест на ветвление в бесконечности при чётном r: если бы бесконечность
#     при r = 4 не ветвилась, Риман–Гурвиц дал бы g = 1, а не 5.
import sys, json, time
from sage.all import *

out = {}

# (а) Кани–Розен против формулы (2) Codex: g = 1 + 2^(r-2) (r-3)
def g_hyp(m):          # род y^2 = произведение m различных линейных множителей
    return (m - 1) // 2
kr = {}
for r in range(2, 13):
    s = sum(binomial(r, m) * g_hyp(m) for m in range(1, r + 1))
    f = 1 + 2 ** (r - 2) * (r - 3)
    kr[r] = (int(s), int(f), s == f)
print("Кани–Розен vs формула (2):", kr)
out['kani_rosen_vs_formula'] = kr
assert all(v[2] for v in kr.values())

# (б) прямой род проективной модели u_i^2 = t - r_i, t = (u_0/z)^2 ... общая модель: u_i^2 = alpha_i u_0^2 + beta_i z^2?
# Берём модель в P^r: u_0^2 = t z^2 не однородна; используем u_i^2 - u_0^2 + r_i z^2 = 0 при f_0 = t (r_0 = 0),
# т.е. корни {0, r_1, ..., r_{r-1}}; при корнях без нуля сдвигаем t.
def curve_genus(roots):
    roots = [QQ(x) for x in roots]
    sh = roots[0]
    roots = [x - sh for x in roots]      # первый корень -> 0
    n = len(roots)
    P = ProjectiveSpace(QQ, n, names=['z'] + ['u%d' % i for i in range(n)])
    z = P.gens()[0]; u = P.gens()[1:]
    eqs = [u[i] ** 2 - u[0] ** 2 + roots[i] * z ** 2 for i in range(1, n)]
    C = Curve(eqs, P)
    t0 = time.time()
    g = C.genus()
    sm = C.is_smooth()
    return int(g), bool(sm), round(time.time() - t0, 1)

tests = {
    'r=2 roots 0,1': [0, 1],
    'r=3 roots 0,1,-1 (конгруэнтная)': [0, 1, -1],
    'r=4 Саллоус a=4420,b=-3360: 0,b,-b,-2a': [0, -3360, 3360, -8840],
    'r=4 слой A, p0=1/3: -1,1,-4/3,4/3': [-1, 1, QQ(-4) / 3, QQ(4) / 3],
    'r=4 малые корни 0,1,-1,-5': [0, 1, -1, -5],
    'r=4 ВЫРОЖДЕННЫЙ p0=-1: -1,1,0,0 (контроль: не 5)': [-1, 1, 0, 0],
}
res = {}
for name, roots in tests.items():
    try:
        res[name] = curve_genus(roots)
    except Exception as e:
        res[name] = 'ошибка: %s' % e
    print(name, '->', res[name]); sys.stdout.flush()
out['direct_genus'] = res

# (в) r=5: род 17 по формуле; прямой счёт (P^5, 4 квадрики) — если успеет
try:
    alarm(600)
    res5 = curve_genus([0, 1, -1, 2, -3])
    cancel_alarm()
    print('r=5 roots 0,1,-1,2,-3 ->', res5)
    out['direct_genus_r5'] = res5
except Exception as e:
    cancel_alarm()
    print('r=5: не досчитано:', e)
    out['direct_genus_r5'] = 'не досчитано'

# (г) локальная картина на бесконечности: sqrt(1 - r s) в Q[[s]] существует (биномиальный ряд),
# значит sqrt(t - r) = s^{-1/2} * (единица), и все радикалы дают одно расширение Q((s))(sqrt s): индекс 2.
Rs = PowerSeriesRing(QQ, 's', default_prec=8)
s = Rs.gen()
chk = {}
for r in [1, -3360, 3360, -8840]:
    q = (1 - r * s).sqrt()
    chk[r] = bool((q ** 2 - (1 - r * s)).is_zero())   # с точностью O(s^8)
print('sqrt(1 - r s) в Q[[s]]:', chk)
out['sqrt_1_minus_rs_in_Q[[s]]'] = chk

json.dump(out, open('/home/kep/magicKube/fable_review_codex/g1_genus_lemma.json', 'w'), indent=1, ensure_ascii=False, default=str)
print('OK')
