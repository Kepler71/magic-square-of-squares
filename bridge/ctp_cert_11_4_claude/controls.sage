# -*- coding: utf-8 -*-
# Шаг 2: КОНТРОЛИ. Симметрия, линейность, нулевая диагональ, образ E(Q)/2E(Q),
#        воспроизведение опубликованного примера Фишера 3.4 (кривая 571a1).
# Ранг E(Q) нигде не используется.
import sys, json, time
load('/home/kep/magicKube/bridge/ctp_cert_11_4_claude/ctp_quartic_snapshot.sage')

D = '/home/kep/magicKube/bridge/ctp_cert_11_4_claude/'
(m, n, s, b, e, c, rootsM, I, J, phis, order, delta_trip) = load(D + 'setup.sobj')
R.<X> = QQ[]
fM = prod(X - r for r in rootsM)
M = EllipticCurve([0, 0, 0, fM.coefficients(sparse=False)[1], fM.coefficients(sparse=False)[0]])
F = FisherCTP(QQ, I, J, verbose=False)
Qcache = load(D + 'quartics.sobj')
OUT = {}

def mkdelta(trip):
    return F.E.from_comps([QQ(trip[i]) for i in order])

def sqfree(a):
    a = QQ(a); num = a.numerator()*a.denominator()
    sgn = -1 if num < 0 else 1; num = abs(num); red = sgn
    for p2, k2 in factor(num):
        if k2 % 2:
            red *= p2
    return red

def mul(t1, t2):
    return tuple(sqfree(QQ(a)*QQ(b2)) for a, b2 in zip(t1, t2))

def quartic(trip):
    key = str(tuple(trip))
    if key in Qcache:
        return [QQ(x) for x in Qcache[key]]
    t0 = time.time()
    g = F.quartic_from_delta(mkdelta(trip), tag=key)
    print("  построена квартика %s за %.1f c" % (key, time.time()-t0)); sys.stdout.flush()
    Qcache[key] = [str(x) for x in g]
    save(Qcache, D + 'quartics.sobj')
    return g

TARGET = (1, 274, 274)
g1 = quartic(TARGET)
print("g1 =", g1)

# ================================================ КОНТРОЛЬ 0: пример Фишера 3.4 (571a1), над Q
print("\n=== КОНТРОЛЬ 0: опубликованный пример Фишера 3.4 (571a1) ===")
Ef = EllipticCurve([0, -1, 1, -929, -10595])
If, Jf = IJ_of_curve(Ef)
print("I,J =", If, Jf, " (в статье 44608, 18842960)")
Ff = FisherCTP(QQ, If, Jf, verbose=False)
print("Delta =", Ff.disc, " (в статье -2338816)")
G1 = [1, 0, -8, 20, -2]; G2 = [2, 0, -4, 4, 7]; G3 = [-1, 0, 34, -12, -47]
try:
    gamf, mf = Ff.gamma1(G1, G2, G3)
    print("m =", mf)
    print("gamma1 =", gamf, " (в статье 4/9*(5x^2-16xz-12z^2) = [20/9,-64/9,-16/3])")
    vf, nf = Ff.pair(G1, G2, G3, reps=2)
    vfr, _ = Ff.pair(G2, G1, G3, reps=2)
    print("<g1,g2> =", vf, "   <g2,g1> =", vfr, "  (в статье: нетривиально)")
    OUT['fisher_example'] = {'I': str(If), 'J': str(Jf), 'gamma1': [str(x) for x in gamf],
                             'pair': int(vf), 'pair_reverse': int(vfr)}
except Exception as ex:
    print("ОШИБКА в контроле 0:", repr(ex))
    OUT['fisher_example'] = {'error': repr(ex)}
sys.stdout.flush()

# ================================================ КОНТРОЛЬ 1: образ E(Q)/2E(Q)
print("\n=== КОНТРОЛЬ 1: классы РАЦИОНАЛЬНЫХ точек E(Q) должны спариваться с g1 в НОЛЬ ===")
e1, e2, e3 = rootsM
pts_classes = []
# точки порядка 2 (даны явно, без всякого поиска и без границы ранга)
pts_classes.append(('2-кручение (e1,0)', ((e1-e2)*(e1-e3), e1-e2, e1-e3)))
pts_classes.append(('2-кручение (e2,0)', (e2-e1, (e2-e1)*(e2-e3), e2-e3)))
pts_classes.append(('2-кручение (e3,0)', (e3-e1, e3-e2, (e3-e1)*(e3-e2))))
# поиск нетривиальных точек прямым перебором x (это НЕ граница ранга: точка — это точка)
fx = fM
found_pts = []
for x0 in range(-3000000, 3000000, 30):
    v = fx(x0)
    if v > 0 and v.is_square():
        found_pts.append(x0)
        if len(found_pts) >= 4:
            break
print("найденные x рациональных точек (перебор):", found_pts)
for x0 in found_pts:
    if x0 in (e1, e2, e3):
        continue
    pts_classes.append(('точка x=%d' % x0, (x0-e1, x0-e2, x0-e3)))

ctrl1 = []
for name, trip in pts_classes:
    tr = tuple(sqfree(t) for t in trip)
    t3 = mul(TARGET, tr)
    print("\n  %s: класс %s" % (name, str(tr)))
    if tr == (1, 1, 1):
        print("    тривиальный класс, пропуск"); continue
    try:
        g2 = quartic(tr); g3 = quartic(t3)
        v, npl = F.pair(g1, g2, g3, reps=2)
        vr, _ = F.pair(g2, g1, g3, reps=2)
        print("    <g1,h> = %s   <h,g1> = %s   (обязано быть 0)" % (v, vr))
        ctrl1.append({'name': name, 'class': [str(t) for t in tr], 'pair': int(v), 'pair_rev': int(vr)})
    except Exception as ex:
        print("    ОШИБКА:", repr(ex))
        ctrl1.append({'name': name, 'class': [str(t) for t in tr], 'error': repr(ex)})
    sys.stdout.flush()
OUT['control_image_EQ'] = ctrl1

# ================================================ КОНТРОЛЬ 2: линейность + симметрия по всей Sel^2
print("\n=== КОНТРОЛЬ 2: <g1,.> на ВСЕЙ группе Sel^2: линейность, симметрия, диагональ ===")
SEL = [(-28770,2,-14385),(-28770,137,-210),(-2055,210,-1918),(-2055,14385,-7),
       (-959,1,-959),(-959,274,-14),(-274,105,-28770),(-274,28770,-105),
       (-105,2,-210),(-105,137,-14385),(-30,210,-7),(-30,14385,-1918),
       (-14,1,-14),(-14,274,-959),(-1,105,-105),(-1,28770,-28770),
       (1,1,1),(1,274,274),(14,105,30),(14,28770,2055),
       (30,2,15),(30,137,4110),(105,210,2),(105,14385,137),
       (274,1,274),(274,274,1),(959,105,2055),(959,28770,30),
       (2055,2,4110),(2055,137,15),(28770,210,137),(28770,14385,2)]
vals = {}
detail = {}
for tr in SEL:
    if tr == (1, 1, 1):
        vals[tr] = 0; continue
    t3 = mul(TARGET, tr)
    try:
        g2 = quartic(tr); g3 = quartic(t3)
        if g2[0] == 0 and g3[0] == 0:
            vals[tr] = 0; detail[str(tr)] = 'trivial'; continue
        v, npl = F.pair(g1, g2, g3, reps=2)
        vr, _ = F.pair(g2, g1, g3, reps=2)
        v13, _ = F.pair(g1, g3, g2, reps=2)     # Fisher 3.2(v): <g1,g2> = <g1,g3>
        vals[tr] = int(v)
        detail[str(tr)] = {'pair': int(v), 'rev': int(vr), 'g1g3': int(v13), 'places': int(npl)}
        flag = ''
        if vr != v: flag += '  !!! НАРУШЕНА СИММЕТРИЯ'
        if v13 != v: flag += '  !!! НАРУШЕНО 3.2(v)'
        print("  %-26s <g1,h>=%d  <h,g1>=%d  <g1,g3>=%d%s" % (str(tr), v, vr, v13, flag))
    except Exception as ex:
        print("  %-26s ОШИБКА %s" % (str(tr), repr(ex)))
        detail[str(tr)] = {'error': repr(ex)}
    sys.stdout.flush()

# линейность: <g1, u*w> = <g1,u> + <g1,w>
print("\n  проверка линейности <g1, u*w> = <g1,u> + <g1,w>:")
bad = 0; ok = 0
for u in SEL:
    for w in SEL:
        uw = mul(u, w)
        if u in vals and w in vals and uw in vals:
            if (vals[u] + vals[w]) % 2 != vals[uw]:
                bad += 1
                if bad <= 5:
                    print("    НАРУШЕНИЕ: u=%s w=%s uw=%s : %d+%d != %d" % (u, w, uw, vals[u], vals[w], vals[uw]))
            else:
                ok += 1
print("    линейность: пройдено %d, нарушений %d" % (ok, bad))
OUT['linearity_ok'] = int(ok); OUT['linearity_bad'] = int(bad)
OUT['pairing_on_Sel'] = {str(k): int(v) for k, v in vals.items()}
OUT['detail'] = detail
print("\n  <g1,g1> (диагональ, вычислена) =", vals.get(TARGET))
OUT['diagonal_g1_g1'] = int(vals.get(TARGET, -1))

json.dump(OUT, open(D + 'controls.json', 'w'), indent=1, ensure_ascii=False)
print("\nЗАПИСАНО controls.json")
