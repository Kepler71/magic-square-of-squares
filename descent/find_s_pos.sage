# поиск кривых Кремоны с ненулевым рангом спаривания (PARI ellrank[2] = s > 0)
import os
lim = int(os.environ.get('CMAX', '600')); found = []
for cond in range(11, lim):
    try: lst = list(cremona_curves([cond]))
    except Exception: continue
    for E in lst:
        try:
            rk = pari(E).ellinit().ellrank()
            s = ZZ(rk[2])
            if s > 0:
                nc = len(pari(E).ellinit().ell2cover())
                found.append((E.cremona_label(), int(s), int(nc)))
                print(E.cremona_label(), 's =', int(s), 'покрытий', nc, flush=True)
        except Exception: pass
    if len(found) >= 8: break
print('найдено', len(found))
