from uniform_check import *
r, s = 204, 247
pats = [('s', 'd', 'u'), ('-s', '-d', '-u'), ('r', 'd', '-u'), ('-r', '-d', 'u')]
primes = [29, 31, 37]
rows = []
for (r2, s2, tag) in [(204, 247, 'база'), (2239, 2874, '≡ mod 33263'), (2662, 2897, '≡ mod 33263'), (43, 60, '≢'), (59, 60, '≢'), (52, 61, '≢'), (61, 67, '≢')]:
    line = []
    for l in primes:
        S, good = run(r2, s2, pats, [l])
        line.append(S)
    rows.append((r2, s2, tag, line))
base = rows[0][3]
for r2, s2, tag, line in rows:
    print(f'{r2}/{s2} ({tag}): ' + '; '.join(f'ℓ={l}: |A|={len(S)}, =база: {S == B}' for l, S, B in zip(primes, line, base)), flush=True)
