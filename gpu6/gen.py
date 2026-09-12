"""Таблица различных D(t) = A(t)^2 - 1 = 4ab(b^2-a^2)/(a^2+b^2)^2, t = a/b, gcd(a,b)=1, 2<=b<=B, 0 < a/b < sqrt2-1
(та же сетка, что rational_motion/scripts/grid.py у Codex), вычеты по простым, маски условия
«1+x+y и 1-x-y — квадраты mod p» (как six_lines_exact/scripts/prepare.py у Codex; 0 считается квадратом).
Экономная по памяти версия (numpy int64: N, D < 2^63 при B <= 8192).
Выход: gen_B.bin (n x k байт вычетов; p = «не фильтровать», если p | знаменатель), gen_B_ab.npy (a, b), gen_B_meta.txt.
WIDE=1 — 30 простых < 128 (иначе 17 простых < 64)."""
import sys, os
from math import isqrt
import numpy as np

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61] + \
         ([67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127] if os.environ.get('WIDE') else [])

def main(B):
    assert B <= 8192
    As, Bs = [], []
    for b in range(2, B + 1):
        amax = isqrt(2*b*b) - b
        a = np.arange(1, amax + 1, dtype=np.int64)
        a = a[np.gcd(a, b) == 1]
        As.append(a); Bs.append(np.full(len(a), b, dtype=np.int64))
    a = np.concatenate(As); b = np.concatenate(Bs); del As, Bs
    N = 4*a*b*(b*b - a*a); D = (a*a + b*b)**2                   # < 2^63 при B <= 8192
    assert N.max() < 2**62 and D.max() < 2**62
    g = np.gcd(N, D); N //= g; D //= g
    order = np.argsort(N.astype(np.float64) / D.astype(np.float64), kind='stable')
    a, b, N, D = a[order], b[order], N[order], D[order]
    # различность: соседние по float — точное сравнение (Python int)
    fl = N.astype(np.float64) / D.astype(np.float64)
    close = np.nonzero(np.diff(fl) <= 1e-12 * np.abs(fl[1:]))[0]
    for i in close:
        # сравнить все элементы окна равных float точно
        assert int(N[i]) * int(D[i+1]) != int(N[i+1]) * int(D[i]), "duplicate D values"
    n = len(a)
    res = np.zeros((n, len(PRIMES)), dtype=np.uint8)
    for l, p in enumerate(PRIMES):
        dm = D % p; nm = N % p
        inv = np.array([0] + [pow(int(x), -1, p) for x in range(1, p)], dtype=np.int64)
        r = (nm * inv[dm]) % p
        r[dm == 0] = p
        res[:, l] = r.astype(np.uint8)
    masks = []
    for p in PRIMES:
        sq = {i*i % p for i in range(p)}
        m = [sum(1 << y for y in range(p) if (1 + x + y) % p in sq and (1 - x - y) % p in sq) | (1 << p) for x in range(p)]
        m += [(1 << (p + 1)) - 1]
        masks.append(m + [0]*(128 - len(m)))
    res.tofile(f'gen_{B}.bin')
    np.save(f'gen_{B}_ab.npy', np.array([a, b], dtype=np.int64))
    with open(f'gen_{B}_meta.txt', 'w') as f:
        f.write(f'{n} {len(PRIMES)}\n' + ' '.join(map(str, PRIMES)) + '\n')
        for m in masks: f.write(' '.join(f'{v & (2**64-1)} {v >> 64}' for v in m) + '\n')   # 128-битные маски
    print(B, n, 'close-float pairs checked exactly:', len(close))

if __name__ == '__main__':
    main(int(sys.argv[1]))
