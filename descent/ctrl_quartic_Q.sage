# Дополнительный контроль над Q: сверка ранга спаривания с PARI ellrank.
# PARI ellrank(E) = [r, R, s, L]; s — ранг спаривания Касселса–Тейта на Sel^2 (PARI реализует ту же
# формулу Фишера, ellrank.c: casselspairing). Квартики базиса Sel^2 берём из PARI ell2cover,
# квартики сумм строим своей конструкцией quartic_from_delta (то есть проверяется и она).
# Запуск: N=40 SEED=1 sage ctrl_quartic_Q.sage
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


def run(E, verbose=False):
    (rlo, rhi, s), gs = covers_of(E)
    n = len(gs)
    if n < 2:
        return ('skip', n, s)
    I, J = IJ_of_curve(E)
    F = FisherCTP(QQ, I, J)
    for g in gs:
        F.check_quartic(g)
    M = matrix(GF(2), n, n)
    for i in range(n):
        for j in range(i + 1, n):
            d3 = F.z_inv(gs[i]) * F.z_inv(gs[j])
            g3 = F.quartic_from_delta(d3, tag=f"{i}_{j}")
            v, npl = F.pair(gs[i], gs[j], g3)
            M[i, j] = v; M[j, i] = v
    gtr = F.quartic_from_delta(F.L(1), tag="triv")
    for i in range(n):
        M[i, i] = F.pair(gs[i], gs[i], gtr)[0]
    return ('case', n, s, M.rank(), bool(M == M.transpose()),
            all(M[i, i] == 0 for i in range(n)), rlo, rhi)


if __name__ == '__main__' or True:
    random.seed(int(os.environ.get('SEED', '1')))
    N = int(os.environ.get('N', '40'))
    rows = []
    bad = 0
    cnt = 0
    for lbl in sorted(set(cremona_curves(srange(1, 600)))) if False else []:
        pass
    E_list = []
    for cond in range(11, 900):
        try:
            for E in cremona_curves([cond]):
                E_list.append(E)
        except Exception:
            pass
    random.shuffle(E_list)
    for E in E_list:
        if cnt >= N:
            break
        try:
            res = run(E)
        except Exception as ex:
            print(f"ERR {E.cremona_label() if hasattr(E,'cremona_label') else E.ainvs()}: "
                  f"{type(ex).__name__}: {str(ex)[:120]}")
            continue
        if res[0] == 'skip':
            continue
        _, n, s, rk, sym, dg0, rlo, rhi = res
        ok = (rk == s)
        if not ok:
            bad += 1
        cnt += 1
        print(f"{E.cremona_label()}: dim(ell2cover) = {n}, PARI s = {s}, мой rank CTP = {rk}, "
              f"{'OK' if ok else 'РАСХОЖДЕНИЕ!!!'}  симм {sym} диаг0 {dg0}  (ранг {rlo}..{rhi})")
        rows.append((E.cremona_label(), n, s, rk, ok, sym, dg0))
    print()
    print(f"проверено {cnt} кривых, расхождений {bad}")
