# Fable, 14.09.2026. CAS-контроль (не часть доказательства): ранги пяти эллиптических факторов якобиана кривой рода 5
# для Саллоуса, корни {0, 3360, -3360, -8840}. PARI ellrank (2-спуск) с лимитом времени на кривую.
import json, itertools, time
from sage.all import *
roots = [0, 3360, -3360, -8840]
out = {}
def quartic_to_weierstrass(rs):
    # y^2 = prod (t - r_i), 4 корня; t = r4 + 1/u => Y^2 = u * prod_{i<4} ((r4 - r_i) u + 1)
    r4 = rs[3]
    Ru = PolynomialRing(QQ, 'u'); u = Ru.gen()
    # (y u^2)^2 = prod_{i<4} ((r4 - r_i) u + 1): кубика c3 u^3 + c2 u^2 + c1 u + c0; X = c3 u, Y = c3 y u^2:
    # Y^2 = X^3 + c2 X^2 + c1 c3 X + c0 c3^2   (исправлено: раньше был лишний множитель u)
    f = prod((r4 - r) * u + 1 for r in rs[:3])
    c = f.coefficients(sparse=False)
    c0, c1, c2, c3 = c[0], c[1], c[2], c[3]
    E = EllipticCurve([0, c2, 0, c1 * c3, c0 * c3 ** 2])
    # контроль: у квартики y^2 = prod(t - r_i) есть точки (r_i, 0) и две на бесконечности; проверим, что 2-кручение E полное
    assert E.torsion_subgroup().order() % 4 == 0
    return E
curves = {}
for T in itertools.combinations(roots, 3):
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    f = prod(x - r for r in T)
    cf = f.coefficients(sparse=False)
    curves['cubic %s' % (T,)] = EllipticCurve([0, cf[2], 0, cf[1], cf[0]])
curves['quartic all four'] = quartic_to_weierstrass(roots)
for name, E in curves.items():
    Em = E.minimal_model()
    t0 = time.time()
    try:
        alarm(400)
        rk = pari(Em).ellrank()
        cancel_alarm()
        res = {'lower': int(rk[0]), 'upper': int(rk[1]), 'cond_digits': len(str(Em.conductor())), 'sec': float(round(time.time() - t0, 1))}
    except Exception as e:
        cancel_alarm()
        res = {'error': str(e)[:80], 'sec': float(round(time.time() - t0, 1))}
    out[name] = res
    print(name, res, flush=True)
json.dump(out, open('/home/kep/magicKube/fable_review_codex/g4_ranks.json', 'w'), indent=1, ensure_ascii=False)
print('OK')
