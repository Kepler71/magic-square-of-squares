# Контроль приёма Фишера на кривых с НЕНУЛЕВЫМ спариванием (PARI s > 0).
import os, time
load('/home/kep/magicKube/descent/ctp_quartic.sage')
load('/home/kep/magicKube/descent/ctp_trick.sage')
def covers_of(E):
    pe = pari(E).ellinit(); rk = pe.ellrank(); C = pe.ell2cover()
    Rx = QQ['x']; gs = []
    for c in C:
        q = Rx(c[0]); gs.append([QQ(q[4 - j]) for j in range(5)])
    return (ZZ(rk[0]), ZZ(rk[1]), ZZ(rk[2])), gs
labels = os.environ.get('LAB', '66b3,102b5,114c3,120a5,130b3,195a7,210c5,210e5').split(',')
NCHK = int(os.environ.get('NCHK', '4'))
tot_pairs = 0; tot_mism = 0; ok_rank = 0; bad_rank = 0
for lab in labels:
    E = EllipticCurve(lab)
    (rlo, rhi, s), gs = covers_of(E); n = len(gs)
    I, J = IJ_of_curve(E); F = FisherCTP(QQ, I, J)
    try:
        for g in gs: F.check_quartic(g)
        M = matrix(GF(2), n, n); t0 = time.time(); npair = 0
        for i in range(n):
            for j in range(i + 1, n):
                # квартика класса-суммы: ищем её среди покрытий PARI (z(g1)z(g2)z(g3) — квадрат в L)
                g3 = None
                for m3 in range(n):
                    if F.E.is_square(F.z_inv(gs[i]) * F.z_inv(gs[j]) * F.z_inv(gs[m3])):
                        g3 = gs[m3]; break
                if g3 is None:
                    print(f"  {lab} ({i},{j}): квартика суммы не найдена среди покрытий"); continue
                old, _ = F.pair(gs[i], gs[j], g3, verbose=False)
                new, used = pair_trick(F, gs[i], gs[j], g3, bound=4, verbose=False, ncheck=NCHK)
                npair += 1; tot_pairs += 1
                pts = [str(u[0]) for u in used]
                if old != new:
                    tot_mism += 1
                    print(f"  {lab} ({i},{j}): РАСХОЖДЕНИЕ старое {old} новое {new}, точки {pts}")
                M[i, j] = new; M[j, i] = new
        good = (M.rank() == s); ok_rank += good; bad_rank += (not good)
        print(f"{lab}: покрытий {n}, PARI s = {s}, мой rank = {M.rank()} {'OK' if good else 'РАСХОЖДЕНИЕ'}, "
              f"пар {npair}, точек на пару {NCHK}  ({time.time()-t0:.0f}s)")
    except Exception as ex:
        print(f"ERR {lab}: {type(ex).__name__}: {str(ex)[:100]}")
print(f"\nИТОГ: пар сверено {tot_pairs}, расхождений старое/новое {tot_mism}; ранг против PARI: OK {ok_rank}, плохо {bad_rank}")
