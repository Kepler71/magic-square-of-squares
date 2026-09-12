# Прогон по всем 16 конфигурациям: 4 типа семейств x 10 базовых пар, решётки Морделла–Вейля,
# LLL + перебор по канонической высоте с потолком + модульное сито + точная проверка выживших.
#   sage scan.sage MODE I0 I1 [K hcap nprimes tmax fams]
#   MODE = sq  — пары (P,Q) с КВАДРАТНЫМ центром (P^2+Q^2 = 2 c0^2), перечисление по (m,n)
#   MODE = gen — все взаимно простые пары P > Q >= 1 (центр вообще говоря не квадрат)
import functools, sys, time, itertools
print = functools.partial(print, flush=True)
argv = sys.argv[1:]; sys.argv = ['x']
load('/home/kep/magicKube/bremner16/engine.sage')

MODE = argv[0]
I0, I1 = int(argv[1]), int(argv[2])
K = int(argv[3]) if len(argv) > 3 else 1500
HCAP = int(argv[4]) if len(argv) > 4 else 1200
NPR = int(argv[5]) if len(argv) > 5 else 60
TMAX = int(argv[6]) if len(argv) > 6 else 120
FAMS = argv[7].split(',') if len(argv) > 7 else ['F1', 'F2', 'F3', 'F4']


def sq_pairs(limit=400):
    """ (P,Q) с P^2+Q^2 = 2*c0^2: P = |m^2-n^2+2mn|, Q = |m^2-n^2-2mn|, c0 = m^2+n^2. """
    out = set()
    for m in range(2, limit):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            P = abs(m^2 - n^2 + 2*m*n); Q = abs(m^2 - n^2 - 2*m*n)
            if P == 0 or Q == 0 or P == Q:
                continue
            g = gcd(P, Q); P, Q = P//g, Q//g
            if P < Q:
                P, Q = Q, P
            if (P^2 + Q^2) % 2 == 0 and ZZ((P^2 + Q^2)//2).is_square():
                out.add((P, Q))
    return sorted(out, key=lambda pq: (pq[0], pq[1]))


def gen_pairs(limit):
    """ Взаимно простые P > Q >= 1 с НЕквадратным центром (квадратный центр идёт режимом sq). """
    out = [(P, Q) for P in range(2, limit) for Q in range(1, P)
           if gcd(P, Q) == 1 and not ((P^2 + Q^2) % 2 == 0 and ZZ((P^2 + Q^2)//2).is_square())]
    return sorted(out, key=lambda pq: (pq[0], pq[1]))


pairs = sq_pairs(400) if MODE == 'sq' else gen_pairs(200)
pairs = pairs[I0:I1]
print(f"# MODE={MODE} пар={len(pairs)} K={K} hcap={HCAP} primes={NPR} tmax={TMAX} fams={FAMS}")

TOTAL = {'bases': 0, 'err': 0, 'nocurve': 0, 'timeout': 0, 'states': 0, 'exact': 0, 'hits': 0}
for (P, Q) in pairs:
    C = (QQ(P)^2 + QQ(Q)^2)/2
    csq = C.is_square()
    print(f"== (P,Q)=({P},{Q}) C={C} {'SQ' if csq else 'non'}")
    for ft in FAMS:
        cells, pair, j, free = family_cells(ft, P, Q)
        ts = degenerate_ts(cells)
        base_cnt = 3 + int(csq)
        for c1, c2 in itertools.combinations(free, 2):
            cands = [(t0, x, y) for t0 in ts for (x, y) in ((c1, c2), (c2, c1))
                     if (cells[x][0] + cells[x][1]*t0^2) > 0 and (cells[x][0] + cells[x][1]*t0^2).is_square()]
            if not cands:
                TOTAL['nocurve'] += 1
                print(f"  {ft} {NAMES[c1]}/{NAMES[c2]}: нет рациональной точки для коники")
                continue
            done = False
            for (t0, x, y) in cands:
                try:
                    alarm(TMAX)
                    line, hits = run_base(cells, free, base_cnt, x, y, t0, K=K, nprimes=NPR, hcap=HCAP)
                    cancel_alarm()
                except (AlarmInterrupt, KeyboardInterrupt):
                    cancel_alarm(); TOTAL['timeout'] += 1
                    print(f"  {ft} {NAMES[x]}/{NAMES[y]}@{t0}: TIMEOUT {TMAX}s"); done = True; break
                except Exception as ex:
                    cancel_alarm()
                    if 'квартик' in str(ex) or 'нет точки' in str(ex) or 'вырожден' in str(ex) or 'singular' in str(ex):
                        continue
                    TOTAL['err'] += 1
                    print(f"  {ft} {NAMES[x]}/{NAMES[y]}@{t0}: ERROR {type(ex).__name__}: {str(ex)[:70]}")
                    done = True; break
                TOTAL['bases'] += 1
                print(f"  {ft} " + line)
                for hh in hits:
                    TOTAL['hits'] += 1
                    print(f"   HIT {ft} (P,Q)=({P},{Q}) t={hh[2]} nsq={hh[0]} pos={hh[1]} cells={hh[3]}")
                done = True; break
            if not done:
                TOTAL['nocurve'] += 1
                print(f"  {ft} {NAMES[c1]}/{NAMES[c2]}: кривая не построена (нет точки на квартике)")
print("# ИТОГ", TOTAL)
