# Независимая проверка утверждения 4 (применимость BST Thm 1.1).
# ЧИТАЛ до написания: NOTE.md §2.4 (доказательство) и чужой скрипт parity_bst_applicability.py.
# Здесь проверка устроена иначе: (1) полный перебор ЛЕММ по вычетам (это доказательство, т.к. вычеты квадратов конечны);
# (2) перебор троек (a,b,c) mod 2^k через numpy до k=9; (3) проверка приведения к бесквадратной модели на символах.
import itertools, numpy as np
from fractions import Fraction

# (1) Леммы. Q16 = квадраты по модулю 16.
Q16 = sorted({x*x % 16 for x in range(16)})
print("квадраты mod 16:", Q16)
# Лемма A: центр a -- квадрат, чётный => a ≡ 0 (mod 4); тогда s1+s2 ≡ 2a ≡ 0 (mod 8) влечёт s1,s2 чётны.
bad = [(a, s1, s2) for a in Q16 if a % 2 == 0 for s1 in Q16 for s2 in Q16
       if (s1 + s2 - 2*a) % 16 == 0 and (s1 % 2 or s2 % 2)]
print("лемма A (a чётно => пара чётна), контрпримеров:", len(bad))
# Лемма B: a нечётно => s1+s2 ≡ 2a (mod 16) влечёт s1,s2 нечётны.
bad = [(a, s1, s2) for a in Q16 if a % 2 == 1 for s1 in Q16 for s2 in Q16
       if (s1 + s2 - 2*a) % 16 == 0 and (s1 % 2 == 0 or s2 % 2 == 0)]
print("лемма B (a нечётно => пара нечётна), контрпримеров:", len(bad))
# Лемма C: нечётный квадрат ≡ 1 (mod 8).
print("лемма C: нечётные квадраты mod 8 =", sorted({x*x % 8 for x in range(1, 8, 2)}))
# Четыре пары через центр покрывают все 8 нецентральных клеток {(i,j)} и в каждой (i,j)+(−i,−j):
pairs = [((1,0),(-1,0)), ((0,1),(0,-1)), ((1,1),(-1,-1)), ((1,-1),(-1,1))]
cover = sorted(p for pr in pairs for p in pr)
print("пары через центр покрывают все 8 клеток:", cover == sorted((i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0)))

# (2) Перебор mod 2^k (векторно).
for k in range(3, 10):
    M = 2**k
    issq = np.zeros(M, dtype=bool); issq[[(x*x) % M for x in range(M)]] = True
    sqs = np.nonzero(issq)[0]
    # центр a обязан быть квадратом: перебираем a по квадратам, b,c по всем вычетам
    a = sqs[:, None, None]; b = np.arange(M)[None, :, None]; c = np.arange(M)[None, None, :]
    ok = np.ones((len(sqs), M, M), dtype=bool)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            ok &= issq[(a + i*b + j*c) % M]
    A, B, C = np.nonzero(ok); A = sqs[A]
    odd = (A % 2 == 1)
    stepbad = odd & (((B % 8) != 0) | ((C % 8) != 0) | (((B + C) % 8) != 0) | (((B - C) % 8) != 0))
    # смешанная чётность: есть и чётная, и нечётная клетка
    par = np.stack([((A + i*B + j*C) % M) % 2 for i in (-1,0,1) for j in (-1,0,1)])
    mixed = (par.min(0) != par.max(0))
    print(f"mod 2^{k}: решений {ok.sum()}, нечётный центр {odd.sum()}, смешанная чётность {mixed.sum()}, "
          f"нечётные с шагом ≢0 mod 8: {stepbad.sum()}")

# (3) Приведение к бесквадратной модели: b = N s^2, N бесквадратно; изоморфизм x -> x/s^2.
# Если v2(b) >= 3 и v2(N) <= 1, то v2(s) >= 1; x = a/s^2 при нечётном a имеет v2 < 0.
def v2(n):
    n = abs(n); v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v
viol = 0
for vb in range(3, 12):
    for vN in (0, 1):
        if (vb - vN) % 2: continue
        vs = (vb - vN)//2
        if vs < 1: viol += 1
print("случаи v2(b)>=3 с v2(s)=0:", viol, "(0 => x=a/s^2 никогда не целое при нечётном a)")
# Проверка единственности бесквадратной модели: y^2=x^3-A x изоморфна y^2=x^3-A' x над Q iff A'/A -- четвёртая степень;
# A=b^2, A'=n^2 => n/b = ±u^2 => n = N (бесквадратные части совпадают). Символьно это очевидно; пример:
b = 8*3*25  # v2=3, N=6, s=10
from math import isqrt
N = b; s = 1
for p in range(2, 100):
    while N % (p*p) == 0: N //= p*p; s *= p
print(f"пример b={b}: N={N}, s={s}; x=a/s^2 при a=1: {Fraction(1, s*s)}")
