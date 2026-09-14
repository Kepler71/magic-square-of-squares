# Fable, 14.09.2026. Мелкие сверки: (1) j(E_{1,2}) — уступка по рецензии Codex (CM нет => не изогенна конгруэнтной);
# (2) корни главной диагонали 3P из рецензии Codex; (3) оценка 314505646^2 < N < 314505647^2 для числителя S(2P);
# (4) сверка параметризации B (lambda-набор) с клетками переписи {1, 1±sz, 1±rz, 1±(s−r)z, 1±(s+r)z}.
import json
from sage.all import *
out = {}
# (1) семейство unify: E_{a,b}: y^2 = x(x+a^2)(x+b^2), a=m^2, b=n^2; для (m,n)=(1,2): y^2 = x(x+1)(x+16); также вариант x(x+1)(x+4)
for name, (A, B) in {'x(x+1)(x+16)': (1, 16), 'x(x+1)(x+4)': (1, 4)}.items():
    E = EllipticCurve([0, A + B, 0, A * B, 0])
    j = E.j_invariant()
    out['j ' + name] = {'j': str(j), 'integral': j.is_integer(), 'has_cm': E.has_cm()}
print({k: v for k, v in out.items()})
# (2) 3P
n = 3360
E = EllipticCurve([0, 0, 0, -n ** 2, 0]); P = E(113 ** 2, 127 * 113 * 97)
h3 = (3 * P)[0]
roots3 = [sqrt(QQ(v)) for v in [h3 + n, h3, h3 - n]]
codex3 = [QQ(5789467770809103487) / 68581433798559363, QQ(4208856041891545073) / 68581433798559363, QQ(1382389703917313183) / 68581433798559363]
out['3P_diag_roots_match_codex'] = bool(sorted(roots3) == sorted(codex3))
# (3)
N = 98913801874105761
out['isqrt_bound_2P'] = bool(314505646 ** 2 < N < 314505647 ** 2) and (isqrt(N) == 314505646)
out['gcd_N_den'] = int(gcd(N, 7751179400836))
# (4) наклон k = r/s, p = s z: lambda*p для lambda in {±1, ±k, ±(1+k), ±(1−k)} = {±s z, ±r z, ±(s+r) z, ±(s−r) z}
r, s, z = var('r s z'); k = r / s; p = s * z
lam = [1, -1, k, -k, 1 + k, -1 - k, 1 - k, -1 + k]
setB = sorted(str((l * p).simplify_full()) for l in lam)
setC = sorted(str(v) for v in [s * z, -s * z, r * z, -r * z, (s + r) * z, -(s + r) * z, (s - r) * z, -(s - r) * z])
out['B_matches_census_cells'] = setB == setC
print({k: v for k, v in out.items() if not k.startswith('j ')})
json.dump(out, open('/home/kep/magicKube/fable_review_codex/g5_misc.json', 'w'), indent=1, ensure_ascii=False, default=str)
print('OK')
