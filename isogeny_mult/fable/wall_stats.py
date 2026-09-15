# Fable 15.09. Статистика «стены» по переписи рангов 11 больших классов (ranks_*.jsonl).
# Критерии на наклон: (N) наивный: есть класс с доказанным рангом r < m;  (F) реальный (метод demjanenko2): есть класс с доказанным r ≤ 2;
# (R1) есть большой класс с r ≤ 1 (старый метод §7 на этих классах).  Доказанность: lo == hi.
import json, os, sys
from collections import Counter
D = '/home/kep/magicKube/isogeny_mult/fable'
rows = {}
for f in sorted(os.listdir(D)):
    if f.startswith('ranks_') and f.endswith('.jsonl'):
        for l in open(f'{D}/{f}'):
            try: d = json.loads(l)
            except Exception: continue
            rows[d['slope']] = d['classes']
def stats(smax, smin=3):
    sel = {sl: c for sl, c in rows.items() if smin <= int(sl.split('/')[1]) <= smax}
    N = F = R1 = 0; walls = []; undec = 0; minr = Counter(); walls_F = []
    for sl, cl in sel.items():
        proved = [(m, lo) for m, T, lo, hi in cl if lo is not None and hi is not None and lo == hi]
        if len(proved) < len(cl): undec += 1
        n = any(lo < m for m, lo in proved); fF = any(lo <= 2 for m, lo in proved); r1 = any(lo <= 1 for m, lo in proved)
        N += n; F += fF; R1 += r1
        if not n: walls.append((sl, [(m, lo, hi) for m, T, lo, hi in cl]))
        if not fF: walls_F.append((sl, [(m, lo, hi) for m, T, lo, hi in cl]))
        if proved: minr[min(lo for m, lo in proved)] += 1
    print(f's ∈ [{smin},{smax}]: наклонов {len(sel)}; наивный критерий (r<m) выполнен: {N}; реальный (r≤2): {F}; r≤1 среди 11 классов: {R1}; '
          f'наклонов с недоказанным рангом хотя бы у одного класса: {undec}')
    print('   распределение min доказанного ранга по 11 классам:', sorted(minr.items()))
    print('   стена (наивная) — наклонов:', len(walls), walls[:10])
    print('   стена (реальная, все 11 классов ранга ≥ 3) — наклонов:', len(walls_F), walls_F[:10])
    return walls, walls_F
if __name__ == '__main__':
    print('всего наклонов в переписи рангов:', len(rows))
    stats(200); stats(500, 201); w, wF = stats(500)
    json.dump(dict(naive_walls=w, real_walls=wF), open(f'{D}/wall_stats.json', 'w'))
