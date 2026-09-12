# Прогон по G-семействам (две полные противоположные пары): 4 автоматических квадрата ВСЕГДА,
# базовая пара даёт 5-й и 6-й, седьмая клетка — «везение» среди оставшихся трёх.
#   sage scan2.sage I0 I1 [K hcap nprimes tmax fams]
import functools, sys, time, itertools
print = functools.partial(print, flush=True)
argv = sys.argv[1:]; sys.argv = ['x']
load('/home/kep/magicKube/bremner16/engine2.sage')

I0, I1 = int(argv[0]), int(argv[1])
K = int(argv[2]) if len(argv) > 2 else 1200
HCAP = int(argv[3]) if len(argv) > 3 else 1000
NPR = int(argv[4]) if len(argv) > 4 else 50
TMAX = int(argv[5]) if len(argv) > 5 else 90
FAMS = argv[6].split(',') if len(argv) > 6 else ['G1', 'G2', 'G3']
if len(argv) > 7:
    RATBOUND = int(argv[7])

M = 60
mn = [(m, n) for m in range(1, M + 1) for n in list(range(-M, 0)) + list(range(1, M + 1))
      if gcd(m, abs(n)) == 1 and m != abs(n)]
mn.sort(key=lambda t: (max(t[0], abs(t[1])), t[0], t[1]))
mn = mn[I0:I1]
print(f"# G-прогон: (m,n) {len(mn)} шт., K={K} hcap={HCAP} primes={NPR} tmax={TMAX} fams={FAMS}")

TOT = {'bases': 0, 'err': 0, 'nocurve': 0, 'timeout': 0, 'hits': 0}
for (m, n) in mn:
    print(f"== (m,n)=({m},{n})")
    for gt in FAMS:
        try:
            cells, auto, free4 = cells_G(gt, m, n)
        except Exception as ex:
            TOT['err'] += 1; print(f"  {gt}: ERROR построения {ex}"); continue
        free = free4 + [4]                      # 4 свободные клетки + центр
        ts = degenerate_ts(cells)
        for c1, c2 in itertools.combinations(free, 2):
            done = False
            for t0 in ts:
                for (x, y) in ((c1, c2), (c2, c1)):
                    v = cells[x](t0)
                    if v <= 0 or not v.is_square():
                        continue
                    try:
                        alarm(TMAX)
                        line, hits = run_base(cells, free, 4, x, y, t0, K=K, nprimes=NPR, hcap=HCAP)
                        cancel_alarm()
                    except (AlarmInterrupt, KeyboardInterrupt):
                        cancel_alarm(); TOT['timeout'] += 1
                        print(f"  {gt} {NAMES[x]}/{NAMES[y]}@{t0}: TIMEOUT"); done = True; break
                    except Exception as ex:
                        cancel_alarm()
                        if 'no point' in str(ex) or 'degenerate' in str(ex) or 'singular' in str(ex):
                            continue
                        TOT['err'] += 1
                        print(f"  {gt} {NAMES[x]}/{NAMES[y]}@{t0}: ERROR {type(ex).__name__}: {str(ex)[:60]}")
                        done = True; break
                    TOT['bases'] += 1
                    print(f"  {gt} " + line)
                    for hh in hits:
                        TOT['hits'] += 1
                        print(f"   HIT {gt} (m,n)=({m},{n}) p={hh[2]} nsq={hh[0]} pos={hh[1]} cells={hh[3]}")
                    done = True; break
                if done:
                    break
            if not done:
                TOT['nocurve'] += 1
                print(f"  {gt} {NAMES[c1]}/{NAMES[c2]}: кривая не построена")
print("# ИТОГ", TOT)
