# Третий, чисто комбинаторный путь (без eclib и без PARI в самом вычислении):
# (1) формула Монски для s(n) = dim Sel_2(E_n) - 2 через матрицу над F_2;
#     формула сверяется с PARI ell2cover для всех бесквадратных n <= 600;
# (2) теорема Таннелла (Waldspurger): L(E_n,1) = 0  <=>  A_n = 2 B_n  (безусловно).
import itertools, math, json, sys

def factor(n):
    f, p = [], 2
    while p*p <= n:
        while n % p == 0:
            f.append(p); n //= p
        p += 1
    if n > 1: f.append(n)
    return f

def leg(a, p):   # символ Лежандра (a/p), p нечётное простое
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def add(s):      # аддитивная запись: 1 -> 0, -1 -> 1
    return 0 if s == 1 else 1

def rank_f2(M):
    M = [row[:] for row in M]; r = 0
    rows, cols = len(M), (len(M[0]) if M else 0)
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(rows):
            if i != r and M[i][c]:
                M[i] = [(x ^ y) for x, y in zip(M[i], M[r])]
        r += 1
    return r

def block(A, B, C, D):
    return [a + b for a, b in zip(A, B)] + [c + d for c, d in zip(C, D)]

def madd(X, Y): return [[(x + y) % 2 for x, y in zip(r, s)] for r, s in zip(X, Y)]
def tr(X): return [list(r) for r in zip(*X)]

def monsky_s(n):
    f = factor(n)
    assert len(set(f)) == len(f), "n must be squarefree"
    even = (f and f[0] == 2)
    ps = [p for p in f if p != 2]
    k = len(ps)
    if k == 0:
        return 0  # n = 1, 2: rank 0, Sel_2 = torsion image (проверяется сравнением)
    A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j: A[i][j] = add(leg(ps[j], ps[i]))
        A[i][i] = sum(A[i][j] for j in range(k) if j != i) % 2
    D = lambda u: [[(add(leg(u, ps[i])) if i == j else 0) for j in range(k)] for i in range(k)]
    D2, Dm2, Dm1 = D(2), D(-2), D(-1)
    if not even:
        M = block(madd(A, D2), D2, D2, madd(A, Dm2))
    else:
        M = block(D2, madd(A, D2), madd(tr(A), D2), Dm1)
    return 2*k - rank_f2(M)

def tunnell(n):
    # n бесквадратное; нечётное: 2x^2+y^2+8z^2=n vs 2x^2+y^2+32z^2=n;
    # чётное: 4x^2+y^2+8z^2=n/2 vs 4x^2+y^2+32z^2=n/2.  Возвращает (A, B).
    m, a = (n, 2) if n % 2 else (n // 2, 4)
    R = int(math.isqrt(m)) + 1
    A = B = 0
    for x in range(-R, R+1):
        for z in range(-R, R+1):
            for (cz, which) in ((8, 'A'), (32, 'B')):
                t = m - a*x*x - cz*z*z
                if t < 0: continue
                y = math.isqrt(t)
                if y*y == t:
                    c = 1 if y == 0 else 2
                    if which == 'A': A += c
                    else: B += c
    return A, B

# сверка формулы Монски с PARI ell2cover
pari = {}
for line in open('/home/kep/magicKube/joint_Eb_Ec/verify_compat/sel2_pari.txt'):
    a, b = line.split(); pari[int(a)] = int(b)
bad = [n for n in pari if monsky_s(n) != pari[n] - 2]
print("Monsky vs PARI ell2cover, squarefree n<=600 and 34,3434,374: mismatches =", bad[:20], "count", len(bad), "of", len(pari))
res = {}
for n, label in [(34, 'b=34 (и c=3400=10^2*34)'), (3434, 'b+c=3434'), (374, 'c-b=3366=3^2*374')]:
    s = monsky_s(n); A, B = tunnell(n)
    res[n] = dict(label=label, n_mod_8=n % 8, monsky_s=s, tunnell_A=A, tunnell_B=B, L1_zero=(A == 2*B))
    print(n, res[n])
json.dump(res, open('/home/kep/magicKube/joint_Eb_Ec/verify_compat/monsky_tunnell.json', 'w'), indent=1)
