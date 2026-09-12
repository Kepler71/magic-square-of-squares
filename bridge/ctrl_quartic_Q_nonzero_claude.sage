# Контроль применимости реализации Фишера НАД Q на случаях с НЕНУЛЕВЫМ спариванием.
# Прежний ctrl_quartic_Q.sage проверял только кривые с s = 0 (и падал на нормировке инвариантов).
# Здесь: (1) чиним нормировку квартик из ell2cover (scale = J*I(g)/(J(g)*I), квадрат);
#        (2) берём ТОЛЬКО кривые, где PARI ellrank даёт s >= 1, то есть CTP заведомо ненулевое;
#        (3) сверяем ранг моей матрицы CTP с s.
# Запуск: NEED=8 sage ctrl_quartic_Q_nonzero_claude.sage
import os, time, traceback
load('/home/kep/magicKube/descent/ctp_quartic.sage')


def covers_of(E):
    pe = pari(E).ellinit()
    rk = pe.ellrank()
    C = pe.ell2cover()
    gs = []
    Rx = QQ['x']
    for c in C:
        q = Rx(c[0])
        gs.append([QQ(q[4 - j]) for j in range(5)])
    return (ZZ(rk[0]), ZZ(rk[1]), ZZ(rk[2])), gs


def run(E):
    (rlo, rhi, s), gs0 = covers_of(E)
    if s < 1:
        return None
    I, J = IJ_of_curve(E)
    F = FisherCTP(QQ, I, J)
    gs = []
    for g in gs0:
        scale = QQ(J * quartic_I(g) / (quartic_J(g) * I))
        assert scale.is_square(), 'scale %s не квадрат' % scale
        g = [scale * c for c in g]
        F.check_quartic(g)
        gs.append(F._make_z_unit(g))
    n = len(gs)
    M = matrix(GF(2), n, n)
    for i in range(n):
        for j in range(i + 1, n):
            d3 = F.z_inv(gs[i]) * F.z_inv(gs[j])
            g3 = F.quartic_from_delta(d3, tag='%d_%d' % (i, j))
            v, _ = F.pair(gs[i], gs[j], g3, reps=2)
            w, _ = F.pair(gs[j], gs[i], g3, reps=2)      # симметрия
            assert v == w, 'несимметрично на %s' % E.cremona_label()
            M[i, j] = v; M[j, i] = v
    gtr = F.quartic_from_delta(F.L(1), tag='triv')
    for i in range(n):
        M[i, i] = F.pair(gs[i], gs[i], gtr, reps=2)[0]
    return (n, int(s), int(M.rank()), bool(M == M.transpose()),
            all(M[i, i] == 0 for i in range(n)), int(rlo), int(rhi))


need = int(os.environ.get('NEED', '8'))
done = 0; bad = 0
for cond in range(11, 3000):
    if done >= need:
        break
    try:
        curves = list(cremona_curves([cond]))
    except Exception:
        continue
    for E in curves:
        if done >= need:
            break
        try:
            res = run(E)
        except Exception as ex:
            continue
        if res is None:
            continue
        n, s, rk, sym, dg0, rlo, rhi = res
        ok = (rk == s)
        bad += 0 if ok else 1
        done += 1
        print('%s: dim ell2cover = %d, PARI s = %d, мой rank CTP = %d, %s  симм %s диаг0 %s (ранг %d..%d)'
              % (E.cremona_label(), n, s, rk, 'OK' if ok else 'РАСХОЖДЕНИЕ!!!', sym, dg0, rlo, rhi), flush=True)
print()
print('кривых с s>=1 проверено %d, расхождений %d' % (done, bad))
