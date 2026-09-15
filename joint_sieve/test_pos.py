# Положительный контроль: наклон 247/825 (квадрат Бремнера–Саллоуса), z = 168/425², шесть квадратных клеток.
# Решето по подмножествам этих шести клеток НЕ должно выбрасывать класс известной точки (и его 2^6 знаковых вариантов).
from jsieve import *
import sys
r, s = 247, 825
z = QQ(168) / QQ(425)**2
cells6 = [s, -s, -r, -(s + r), s - r, -(s - r)]
roots = {}
for lam in cells6:
    v = 1 + lam * z; assert v > 0 and v.is_square(), lam
    roots[lam] = v.sqrt()
print('корни', {k: str(v) for k, v in roots.items()})
subsets = [list(T) for k in (3, 4) for T in itertools.combinations(cells6, k)]
t0 = time.time()
facs = []
for T in subsets:
    try:
        f = Factor(T, verbose=True)
    except Exception as e:
        print('пропуск', T, e); continue
    facs.append(f)
print('множителей', len(facs), 'время', time.time() - t0, flush=True)
# выбираем доказанные ранги, сначала малые
sel = sorted([i for i, f in enumerate(facs) if f.rank_proved], key=lambda i: (facs[i].rank, i))
nf = int(sys.argv[1]) if len(sys.argv) > 1 else 4
sel = sel[:nf]
F = [facs[i] for i in sel]
for f in F: print('используем', f.name, 'ранг', f.rank, 'кручение', f.tinv, 'z(P_known)=', f.z_of(f.point(z, prod(roots[t] for t in f.T))))
JS = JointSieve(r, s, cells6, F)
known = JS.sign_classes(z, roots)
print('классы известной точки (базовый):', known[0])
known0 = JS.sign_classes(QQ(0), {lam: QQ(1) for lam in cells6})   # вырожденная точка z=0 — тоже решение ослабленной задачи
known_all = known + known0
N0 = {i: lcm(f.tinv) for i, f in enumerate(F)}
primes = JS.good_primes(3, int(sys.argv[2]) if len(sys.argv) > 2 else 300, range(len(F)))
print('простых', len(primes))
t0 = time.time()
JS.verbose = False
cand, idx, cols, N, log = JS.run(list(range(len(F))), N0, primes, lifts=[(i, 2) for i in range(len(F))] + [(i, 3) for i in range(len(F))], known=known)
for st in log: print(f"  {st['step']}: выживших {st['survivors']} (N={[st['N'][i] for i in range(len(F))]})")
print('итог: выживших', cand.shape[0], 'N=', N, 'время', time.time() - t0)
S = set(map(tuple, cand.tolist()))
kn = set(JS.class_tuple(cl, idx, N) for _, cl in known)
kn0 = set(JS.class_tuple(cl, idx, N) for _, cl in known0)
print('классов известной точки (различных):', len(kn), 'все в выживших:', kn <= S)
print('классов точки z=0 (различных):', len(kn0), 'все в выживших:', kn0 <= S)
print('выжившие = известные ∪ z=0:', S == (kn | kn0), '; лишних:', len(S - kn - kn0))
