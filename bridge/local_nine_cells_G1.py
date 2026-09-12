# -*- coding: utf-8 -*-
# Локальная разрешимость систем «k клеток G1 — квадраты» во всех Q_p.
# Строгая процедура: рекурсия по дискам t = r (mod p^k) + явный Гензель.
# Ветвь закрывается ТОЛЬКО когда значение определено (v_p(f(r)) < k) и не квадрат.
# Если после KMAX уровней ветвь не решена — возвращается UNDECIDED (и это НЕ доказательство).
#
# Запуск:  python3 /home/kep/magicKube/bridge/local_nine_cells_G1.py
from fractions import Fraction as Fr
from math import gcd
import sys

def cellpolys(m, n):
    """пять неавтоматических клеток как однородные квадратичные формы от (a,b), t = a/b.
       Возвращает список функций f(a,b) -> целое (всё умножено на 2, чтобы уйти от 1/2;
       множитель 2 меняет класс квадратов, поэтому вместо этого держим Fraction)."""
    def F0(a, b): return Fr(m*m*b*b + n*n*a*a)
    def F4(a, b): return Fr((m*m + n*n)*(a*a + b*b), 2)
    def F8(a, b): return Fr(n*n*b*b + m*m*a*a)
    def C2(a, b): return F4(a, b) - 2*m*n*a*b
    def C6(a, b): return F4(a, b) + 2*m*n*a*b
    return [F0, F4, F8, C2, C6]

def vp(x, p):
    if x == 0: return None            # бесконечная оценка
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0: num //= p; v += 1
    while den % p == 0: den //= p; v -= 1
    return v

def is_sq(x, p):
    """x in Q, x != 0: квадрат ли в Q_p"""
    v = vp(x, p)
    if v % 2: return False
    y = x / Fr(p)**v
    if p == 2:
        return (y.numerator * y.denominator) % 8 == 1
    u = (y.numerator * pow(y.denominator, p - 2, p)) % p
    return pow(u, (p - 1)//2, p) == 1

# порог определённости: значение f(r) определяет класс квадратов на диске
# t = r + p^k Z_p, если v_p(f(r)) < k - e, где e = v_p(2) (=0 для нечётных p, =1 для p=2)
# и дополнительно для p=2 нужен запас 3 (класс определяется mod 8).
def decided(val, k, p):
    if val is None: return False
    if p == 2: return val + 3 <= k
    return val < k

def solve_disk(fs, r, k, p, KMAX, budget):
    """есть ли t в диске r + p^k Z_p со всеми f_i квадратами.
       True / False / None(undecided)"""
    if budget[0] <= 0: return None
    budget[0] -= 1
    vals = [f(r, 1) for f in fs]
    undec = False
    for x in vals:
        v = vp(x, p)
        if decided(v, k, p):
            if not is_sq(x, p): return False
        else:
            undec = True
    if not undec:
        return True
    if k >= KMAX: return None
    res = False
    for z in range(p):
        out = solve_disk(fs, r + z * p**k, k + 1, p, KMAX, budget)
        if out is True: return True
        if out is None: res = None
    return res

def locally_solvable(m, n, p, KMAX=7, BUDGET=400000):
    """P^1(Q_p) = {t in Z_p} U {t = 1/u, u in pZ_p}"""
    fs = cellpolys(m, n)
    budget = [BUDGET]
    out1 = solve_disk(fs, 0, 0, p, KMAX, budget)
    if out1 is True: return True
    # окрестность бесконечности: t = 1/u, u in p Z_p  ->  подставляем (a,b) = (1, u)
    gs = [(lambda f: (lambda u, _one: f(1, u)))(f) for f in fs]
    budget2 = [BUDGET]
    out2 = solve_disk(gs, 0, 1, p, KMAX, budget2)   # u in pZ_p  <=> диск 0 + p^1 Z_p
    if out2 is True: return True
    if out1 is None or out2 is None: return None
    return False

def subsystem(m, n, idx):
    fs = cellpolys(m, n)
    return [fs[i] for i in idx]

def loc_solv_sub(m, n, p, idx, KMAX=7, BUDGET=400000):
    fs = subsystem(m, n, idx)
    budget = [BUDGET]
    o1 = solve_disk(fs, 0, 0, p, KMAX, budget)
    if o1 is True: return True
    gs = [(lambda f: (lambda u, _one: f(1, u)))(f) for f in fs]
    b2 = [BUDGET]
    o2 = solve_disk(gs, 0, 1, p, KMAX, b2)
    if o2 is True: return True
    if o1 is None or o2 is None: return None
    return False

PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]

if __name__ == "__main__":
    THREE = [0,1,2]        # F0, F4, F8 — семь квадратов (кривая C рода 5)
    FIVE  = [0,1,2,3,4]    # + красные c2, c6 — девять квадратов

    print("="*100)
    print("A. Самотест процедуры: пары с ЯВНОЙ глобальной точкой должны быть всюду локально разрешимы")
    print("="*100)
    # (15,8): m^2+n^2 = 289 = 17^2, t = 1 даёт все девять квадратов (вырожденно)
    for (m,n) in [(15,8),(4,3),(7,1)]:
        bad3 = [p for p in PRIMES if loc_solv_sub(m,n,p,THREE) is not True]
        bad5 = [p for p in PRIMES if loc_solv_sub(m,n,p,FIVE)  is not True]
        print(f"  (m,n)=({m},{n}): m^2+n^2 = {m*m+n*n} (квадрат: {int((m*m+n*n)**0.5)**2 == m*m+n*n});"
              f" не-True для 3 клеток: {bad3}; для 5 клеток: {bad5}")

    print()
    print("="*100)
    print("B. Четыре оставшиеся пары проекта — обе системы, все p <= 113")
    print("="*100)
    for (a,b) in [(13,8),(15,8),(16,5),(19,16)]:
        m,n = (b,a) if a % 2 == 0 else (a,b)
        r3 = {p: loc_solv_sub(m,n,p,THREE) for p in PRIMES}
        r5 = {p: loc_solv_sub(m,n,p,FIVE)  for p in PRIMES}
        bad3 = [p for p,v in r3.items() if v is False]; und3 = [p for p,v in r3.items() if v is None]
        bad5 = [p for p,v in r5.items() if v is False]; und5 = [p for p,v in r5.items() if v is None]
        print(f"  ({a},{b}):")
        print(f"     СЕМЬ квадратов  (F0,F4,F8):   препятствие в p = {bad3 or 'нет'};  не решено: {und3 or 'нет'}")
        print(f"     ДЕВЯТЬ квадратов (+c2,c6):    препятствие в p = {bad5 or 'нет'};  не решено: {und5 or 'нет'}")

    print()
    print("="*100)
    print("C. Все 127 пар: сколько закрывается локально на девяти клетках, сколько на семи")
    print("="*100)
    pairs = [(a,b) for a in range(2,21) for b in range(1,a) if gcd(a,b)==1]
    kill3, kill5, und = [], [], []
    for (a,b) in pairs:
        m,n = (b,a) if a % 2 == 0 else (a,b)
        v3 = None; v5 = None
        for p in PRIMES:
            r = loc_solv_sub(m,n,p,THREE)
            if r is False: v3 = p; break
            if r is None: und.append((a,b,p,'3'))
        for p in PRIMES:
            r = loc_solv_sub(m,n,p,FIVE)
            if r is False: v5 = p; break
            if r is None: und.append((a,b,p,'5'))
        if v3: kill3.append((a,b,v3))
        if v5: kill5.append((a,b,v5))
    print(f"  локально закрыто на СЕМИ клетках:  {len(kill3)} из 127")
    print(f"  локально закрыто на ДЕВЯТИ клетках: {len(kill5)} из 127")
    print(f"  пар, НЕ закрытых локально на девяти: {[x for x in pairs if x not in [(a,b) for a,b,_ in kill5]]}")
    print(f"  неразрешённых узлов (UNDECIDED): {len(und)}")
    if und: print("   примеры:", und[:10])
