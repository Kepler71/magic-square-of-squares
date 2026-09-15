# Положительный контроль метода Демьяненко–Манина: симметричные 4-клеточные наборы с известной невырожденной точкой.
from demjanenko import *
pts = []
for q in range(1, 40):
    for p in range(1, 2 * q):
        if gcd(p, q) != 1: continue
        w2 = 2 * q * q - p * p
        if w2 > 0 and ZZ(w2).is_square(): pts.append(QQ(p) / q)
pts = sorted(set(pts))
print('точек окружности', len(pts))
found = 0; tested = 0
for u1, u2 in itertools.combinations(pts, 2):
    a0, b0 = u1**2 - 1, u2**2 - 1
    if a0 == 0 or b0 == 0 or abs(a0) == abs(b0): continue
    D = lcm(a0.denominator(), b0.denominator()); a, b = ZZ(a0 * D), ZZ(b0 * D); z = QQ(1) / D
    g = gcd(a, b); a, b, z = a // g, b // g, z * g
    S = [a, -a, b, -b]; assert all_square(z, S)
    for T in itertools.combinations(sorted(S), 3):
        if sorted(-x for x in T) == list(T): continue
        try: f = Factor(T)
        except Exception: continue
        if not (f.rank_proved and f.rank == 1): continue
        tested += 1
        res = run(S, T, verbose=False)
        ok = str(z) in res['nondeg_solutions']
        print(f'S={S}, z={z}, T={list(T)}: ĥ(P₀)={res["hP0"]:.2f}, M0={res["M0"]}, кандидатов {res["ncand"]}, невырожденные решения {res["nondeg_solutions"]}, известная точка найдена: {ok}', flush=True)
        assert ok, 'ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ ПРОВАЛЕН'
        found += 1
    if found >= 6: break
print('контролей пройдено', found)
