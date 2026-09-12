"""Этап 2: спаривание Касселса-Тейта для G1 (19,5).  Своя реализация формулы Фишера, Thm 3.1.
Ранг НЕ вычисляется.  ell2cover используется ТОЛЬКО как источник кандидатов g2;
всюду-локальная разрешимость всех трёх квартик проверяется здесь заново и ещё раз в верификаторе.
"""
from sage.all import *
import json, time, random
from pathlib import Path

exec(open('/home/kep/magicKube/bridge/ctp_19_5_claude/discover_19_5.py').read().split('t0 = time.time()')[0])

OUT = Path('/home/kep/magicKube/bridge/ctp_19_5_claude')
random.seed(20260912); set_random_seed(20260912)

g1 = quartic_from_delta(DELTA, 'target')
print('g1 =', g1, flush=True)

# ------------------------------------------------ локальная арифметика (своя, точная)
def val_unit(q, p):
    q = QQ(q); assert q != 0
    num, den = ZZ(q.numerator()), ZZ(q.denominator())
    k = 0
    while num % p == 0: num //= p; k += 1
    while den % p == 0: den //= p; k -= 1
    return k, QQ(num)/QQ(den)

def is_loc_square(q, pl):
    q = QQ(q)
    if q == 0: return False
    if pl == 'real': return q > 0
    p = ZZ(pl); k, uu = val_unit(q, p)
    if k % 2: return False
    r = (uu.numerator()*ZZ(uu.denominator()).inverse_mod(8 if p == 2 else p))
    if p == 2: return r % 8 == 1
    return kronecker(r % p, p) == 1

def my_hilbert(a, bb, pl):
    """символ Гильберта (a,b)_v, своя реализация"""
    a, bb = QQ(a), QQ(bb); assert a != 0 and bb != 0
    if pl == 'real':
        return -1 if (a < 0 and bb < 0) else 1
    p = ZZ(pl)
    va, ua = val_unit(a, p); vb, ub = val_unit(bb, p)
    if p == 2:
        ra = ua.numerator()*ZZ(ua.denominator()).inverse_mod(8) % 8
        rb = ub.numerator()*ZZ(ub.denominator()).inverse_mod(8) % 8
        e = ((ra-1)//2)*((rb-1)//2) + va*((rb*rb-1)//8) + vb*((ra*ra-1)//8)
        return -1 if e % 2 else 1
    ra = ua.numerator()*ZZ(ua.denominator()).inverse_mod(p) % p
    rb = ub.numerator()*ZZ(ub.denominator()).inverse_mod(p) % p
    res = 1
    if (va*vb) % 2 and (p % 4) == 3: res = -res
    if vb % 2 and kronecker(ra, p) == -1: res = -res
    if va % 2 and kronecker(rb, p) == -1: res = -res
    return res

# сверка своего символа Гильберта с реализацией Sage/PARI
for _ in range(400):
    aa = QQ(random.randint(-200, 200))/QQ(random.randint(1, 30))
    bb = QQ(random.randint(-200, 200))/QQ(random.randint(1, 30))
    if aa == 0 or bb == 0: continue
    for pl in ['real', 2, 3, 5, 7, 11, 13, 193]:
        mine = my_hilbert(aa, bb, pl)
        sag = hilbert_symbol(aa, bb, -1 if pl == 'real' else pl)
        assert mine == sag, (aa, bb, pl, mine, sag)
print('OK: свой символ Гильберта совпал с Sage на 400 случайных парах x 8 мест', flush=True)

def supp(q):
    q = QQ(q)
    if q == 0: return set()
    return set(p for p, _ in factor(q.numerator())) | set(p for p, _ in factor(q.denominator()))

# ------------------------------------------------ поиск локальной точки
def local_points(g, pl, need=1, avoid=None, limit=4000):
    """рациональные (x:z) с g(x,z) ненулевым квадратом в Q_v, и avoid(x,z)!=0"""
    res = []
    p = None if pl == 'real' else ZZ(pl)
    cands = [(ZZ(1), ZZ(0)), (ZZ(0), ZZ(1))]
    for i in range(-60, 61):
        cands.append((ZZ(i), ZZ(1)))
    for j in range(-60, 61):
        cands.append((ZZ(1), ZZ(j)))
    if p is not None:
        for k in range(1, 7):
            for i in range(-30, 31):
                cands.append((ZZ(i), p**k)); cands.append((ZZ(i)*p**k, ZZ(1)))
        for i in range(0, min(int(p), 400)):
            cands.append((ZZ(i), ZZ(1)))
    seen = set(); cnt = 0
    for (x, z) in cands:
        if (x, z) == (0, 0): continue
        key = (x, z)
        if key in seen: continue
        seen.add(key); cnt += 1
        if cnt > limit: break
        v = ev(g, QQ(x), QQ(z))
        if v == 0 or not is_loc_square(v, pl): continue
        if avoid is not None and ev(avoid, QQ(x), QQ(z)) == 0: continue
        res.append((QQ(x), QQ(z)))
        if len(res) >= need: return res
    # случайный добор
    for _ in range(60000):
        if p is None:
            x = QQ(random.randint(-500, 500)); z = QQ(random.randint(-500, 500))
        else:
            e = random.randint(-5, 5)
            x = QQ(random.randint(-2000, 2000))*QQ(p)**e
            z = QQ(random.randint(-2000, 2000))
        if x == 0 and z == 0: continue
        v = ev(g, x, z)
        if v == 0 or not is_loc_square(v, pl): continue
        if avoid is not None and ev(avoid, x, z) == 0: continue
        if (x, z) in res: continue
        res.append((x, z))
        if len(res) >= need: return res
    return res

# ------------------------------------------------ формула Фишера
def H_forms(g):
    """список из 3 компонент H (по одной на каждый phi_i), каждая = [H0,H1,H2]"""
    h = hessian(g); out = []
    for ph in phis:
        G = [(4*ph*g[j] + h[j])/3 for j in range(5)]
        H = [G[0], G[1]/2, G[2]/6 + QQ(2)/9*(I - ph*ph)]
        Hsq = [sum(H[i2]*H[j2] for i2 in range(3) for j2 in range(3) if i2+j2 == k) for k in range(5)]
        assert Hsq == [G[0]*v for v in G], 'тождество G(1,0)G=H^2 нарушено'
        out.append(H)
    return out

def gamma_form(gA, gB, gC, signs=(1, 1, 1)):
    zA, zB, zC = z_inv(gA), z_inv(gB), z_inv(gC)
    mm = []
    for i in range(3):
        pr = zA[i]*zB[i]*zC[i]
        assert QQ(pr).is_square(), ('z1z2z3 не квадрат в компоненте', i, pr)
        mm.append(signs[i]*QQ(pr).sqrt())
    Hs = H_forms(gA)
    gam = []
    for j in range(3):
        gam.append(sum((mm[i]/zA[i])*Hs[i][j]/prod(phis[i]-phis[k] for k in range(3) if k != i)
                       for i in range(3)))
    return gam, mm, (zA, zB, zC)

def place_set(gA, gam, a):
    S = set([2, 3, 5, 7]) | supp(DISC) | supp(a)
    for c in list(gA) + list(gam):
        S |= supp(c)          # СВЕРХ-множество к Remark 3.3 (консервативно)
    return sorted(S)

def pair_value(gA, gB, gC, signs=(1, 1, 1), reps=2, verbose=False):
    gam, mm, zs = gamma_form(gA, gB, gC, signs)
    a = gA if False else gB[0]
    assert a != 0
    S = place_set(gA, gam, a)
    runs = []
    for r in range(reps):
        row = []
        for pl in [str(p) for p in S] + ['real']:
            pts = local_points(gA, pl, need=r+1, avoid=gam)
            assert len(pts) >= r+1, ('нет локальной точки', pl, len(pts))
            x, z = pts[r]
            q = ev(gA, x, z); gv = ev(gam, x, z)
            hb = my_hilbert(a, gv, pl)
            row.append({'place': pl, 'x': str(x), 'z': str(z), 'g': str(q),
                        'gamma': str(gv), 'hilbert': int(hb)})
        runs.append(row)
        if verbose:
            print('  run', r, 'нетривиальные:', [v['place'] for v in row if v['hilbert'] == -1], flush=True)
    vals = [prod(v['hilbert'] for v in row) for row in runs]
    assert len(set(vals)) == 1, ('разные наборы локальных точек дали разное!', vals)
    return vals[0], gam, mm, zs, S, runs

# ------------------------------------------------ кандидаты g2 из ell2cover
Rp = PolynomialRing(QQ, 'x'); xp = Rp.gen()
cover = pari(M).ell2cover()
print('ell2cover: ', len(cover), 'квартик', flush=True)

found = None
log = []
for i, c in enumerate(cover):
    t1 = time.time()
    try:
        q = Rp(c[0])
        g2 = [QQ(q[4-j]) for j in range(5)]
        Ig, Jg = quartI(g2), quartJ(g2)
        lam2 = QQ(Jg*I)/QQ(J*Ig)
        assert QQ(lam2).is_square(), ('масштаб не квадрат', lam2)
        g2 = [x/lam2 for x in g2]
        assert quartI(g2) == I and quartJ(g2) == J
        g2 = make_z_unit(g2)
        z1 = z_inv(g1); z2 = z_inv(g2)
        d3 = [sqfree(z1[k]*z2[k]) for k in range(3)]
        g3 = quartic_from_delta(d3, 'sum%d' % i)
        val, gam, mm, zs, S, runs = pair_value(g1, g2, g3, verbose=True)
        row = {'i': i, 'lam2': str(lam2), 'g2': [str(x) for x in g2], 'g3': [str(x) for x in g3],
               'z2class': [str(sqfree(v)) for v in z2], 'pairing': int(val),
               'places': len(S)+1, 'seconds': time.time()-t1}
        print('PAIR', row['i'], 'class(z2)=', row['z2class'], '-> <g1,g2> =', val,
              '(%d мест, %.1f c)' % (len(S)+1, time.time()-t1), flush=True)
        log.append(row)
        if val == -1:
            found = (i, g2, g3, gam, mm, zs, S, runs, lam2)
            break
    except Exception as ex:
        print('ERROR', i, repr(ex), flush=True)
        log.append({'i': i, 'error': repr(ex), 'seconds': time.time()-t1})

(OUT/'pairing_scan.json').write_text(json.dumps(log, indent=2)+'\n')
if found is None:
    print('НЕ НАЙДЕНО ненулевого спаривания среди', len(cover), 'кандидатов', flush=True)
    raise SystemExit(1)

i, g2, g3, gam, mm, zs, S, runs, lam2 = found
print('НАЙДЕНО: i =', i)
print('g2 =', g2)
print('g3 =', g3)
print('gamma =', gam)
(OUT/'found.json').write_text(json.dumps({
    'i': i, 'g1': [str(x) for x in g1], 'g2': [str(x) for x in g2], 'g3': [str(x) for x in g3],
    'gamma': [str(x) for x in gam], 'm': [str(x) for x in mm],
    'places': [str(p) for p in S], 'lam2': str(lam2)}, indent=2)+'\n')
print('DONE STAGE2', flush=True)
