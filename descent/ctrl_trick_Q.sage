# Контроль приёма Фишера над Q (письмо 22.09.2026): для каждой пары покрытий сверяем
# старую формулу (значение в точке (1:0)) и новую (точка выбрана приёмом; ncheck разных точек).
# Дополнительно: ранг матрицы спаривания против PARI ellrank (s).
import os, time
load('/home/kep/magicKube/descent/ctp_quartic.sage')
load('/home/kep/magicKube/descent/ctp_trick.sage')

def covers_of(E):
    pe = pari(E).ellinit(); rk = pe.ellrank(); C = pe.ell2cover()
    Rx = QQ['x']; gs = []
    for c in C:
        q = Rx(c[0]); gs.append([QQ(q[4 - j]) for j in range(5)])
    return (ZZ(rk[0]), ZZ(rk[1]), ZZ(rk[2])), gs

random.seed(int(os.environ.get('SEED', '1')))
N = int(os.environ.get('N', '12')); NCHK = int(os.environ.get('NCHK', '3'))
E_list = []
for cond in range(11, 400):
    try:
        E_list += list(cremona_curves([cond]))
    except Exception:
        pass
random.shuffle(E_list)
cnt = 0; pairs = 0; mism = 0; rk_ok = 0; rk_bad = 0
for E in E_list:
    if cnt >= N: break
    try:
        (rlo, rhi, s), gs = covers_of(E)
    except Exception:
        continue
    n = len(gs)
    if n < 2: continue
    I, J = IJ_of_curve(E); F = FisherCTP(QQ, I, J)
    try:
        for g in gs: F.check_quartic(g)
        M = matrix(GF(2), n, n); t0 = time.time(); loc = 0
        for i in range(n):
            for j in range(i + 1, n):
                d3 = F.z_inv(gs[i]) * F.z_inv(gs[j])
                g3 = F.quartic_from_delta(d3, tag=f"{i}_{j}")
                old, _ = F.pair(gs[i], gs[j], g3, verbose=False)
                new, used = pair_trick(F, gs[i], gs[j], g3, bound=4, verbose=False, ncheck=NCHK)
                pairs += 1; loc += 1
                if old != new:
                    mism += 1; print(f"  {E.cremona_label()} ({i},{j}): старое {old} новое {new}")
                M[i, j] = new; M[j, i] = new
    except Exception as ex:
        print(f"ERR {E.cremona_label()}: {type(ex).__name__}: {str(ex)[:90]}"); continue
    good = (M.rank() == s)
    rk_ok += good; rk_bad += (not good); cnt += 1
    print(f"{E.cremona_label()}: покрытий {n}, PARI s = {s}, мой rank = {M.rank()} "
          f"{'OK' if good else 'РАСХОЖДЕНИЕ'}, пар {loc} ({time.time()-t0:.0f}s)")
print(f"\nИТОГ: кривых {cnt}, пар сверено {pairs}, расхождений старое/новое {mism}; "
      f"ранг против PARI: совпало {rk_ok}, разошлось {rk_bad}")
