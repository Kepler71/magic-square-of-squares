#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ЗАДАЧА: проверить ВЫЧИСЛЕНИЕМ, есть ли у (15,8) дешёвое локальное препятствие
        того же типа, что закрыло (13,8) [p=17] и (16,5) [p=41].

Тест (необходимое условие для полной девятки G1):
  t = a/b, однородные формы (получаются из клеток умножением на b^2 — квадрат):
     F0 = m^2 b^2 + n^2 a^2
     F4 = s (a^2 + b^2),      s = (m^2+n^2)/2
     F8 = n^2 b^2 + m^2 a^2
     L  = F4 - 2 m n a b
     U  = F4 + 2 m n a b
  Для рациональной точки нормируем (a,b) p-целыми, хотя бы одно — единица.
  Тогда все пять значений — квадраты в Q_p и p-целые => редуцируются в квадраты F_p (нуль разрешён).
  Достаточно перебрать все p+1 точек P^1(F_p).

Отрицательный результат («препятствия нет») здесь тоже считается, а не предполагается.
"""
from fractions import Fraction

def primes_upto(N):
    sieve = [True]*(N+1); sieve[0]=sieve[1]=False
    for i in range(2,int(N**0.5)+1):
        if sieve[i]:
            for j in range(i*i,N+1,i): sieve[j]=False
    return [i for i,v in enumerate(sieve) if v]

def survivors_mod_p(m, n, p):
    """Все (a:b) in P^1(F_p), где все пять форм — квадраты в F_p (нуль разрешён).
       Возвращает список меток: 'inf' или значение t = a/b в F_p."""
    if p == 2:
        raise ValueError("p=2 обрабатывается отдельно, 2-адически")
    sq = set((i*i) % p for i in range(p))
    inv2 = pow(2, p-2, p)
    s = ((m*m + n*n) % p) * inv2 % p
    m2, n2, mn2 = m*m % p, n*n % p, (2*m*n) % p
    out = []
    pts = [(t, 1) for t in range(p)] + [(1, 0)]
    for a, b in pts:
        a2, b2, ab = a*a % p, b*b % p, a*b % p
        F0 = (m2*b2 + n2*a2) % p
        F4 = s*(a2 + b2) % p
        F8 = (n2*b2 + m2*a2) % p
        L = (F4 - mn2*ab) % p
        U = (F4 + mn2*ab) % p
        if F0 in sq and F4 in sq and F8 in sq and L in sq and U in sq:
            out.append('inf' if b == 0 else a)
    return out

def two_adic_survivors(m, n, N=12):
    """Точный 2-адический тест: существует ли (a,b) в Z_2^2, не оба чётные,
       с пятью значениями — квадратами в Q_2.
       Работаем с точными рациональными значениями при a,b mod 2^N и решаем
       вопрос квадратности в Q_2 только когда оценка и единица mod 8 определены.
       Возвращает список пар (a mod 2^N, b mod 2^N), для которых нашёлся ЛИФТ,
       для которого все пять значений — квадраты в Q_2 (достаточное условие Гензеля)."""
    M = 1 << N
    def is_sq_Q2(fr):
        if fr == 0:
            return True          # нуль разрешён (может быть редукцией)
        num, den = fr.numerator, fr.denominator
        v = 0
        while num % 2 == 0: num //= 2; v += 1
        while den % 2 == 0: den //= 2; v -= 1
        if v % 2: return False
        u = Fraction(num, den)
        # u — 2-адическая единица; квадрат в Z_2^* <=> u = 1 mod 8
        return (u.numerator * pow(u.denominator, -1, 8)) % 8 == 1
    s = Fraction(m*m + n*n, 2)
    good = []
    for a in range(M):
        for b in range(M):
            if a % 2 == 0 and b % 2 == 0: continue
            F0 = Fraction(m*m*b*b + n*n*a*a)
            F4 = s*(a*a + b*b)
            F8 = Fraction(n*n*b*b + m*m*a*a)
            L = F4 - 2*m*n*a*b
            U = F4 + 2*m*n*a*b
            if all(is_sq_Q2(x) for x in (F0, F4, F8, L, U)):
                good.append((a, b))
    return good

def degenerate_ts(m, n):
    """Рациональные t, где какая-то клетка нулевая или две клетки совпадают —
       для (15,8) список из REVIEW_FOR_CLAUDE_GROK §4."""
    if (m, n) == (15, 8):
        return [Fraction(-23,7), Fraction(-15,8), Fraction(-1), Fraction(-8,15),
                Fraction(-7,23), Fraction(0), Fraction(7,23), Fraction(8,15),
                Fraction(1), Fraction(15,8), Fraction(23,7)]
    return []

if __name__ == '__main__':
    import sys, math
    P = primes_upto(200)
    m, n = 15, 8
    print("="*78)
    print("1. (15,8): ПОЛНЫЙ ПЕРЕБОР P^1(F_p) ДЛЯ ВСЕХ НЕЧЁТНЫХ p < 200")
    print("="*78)
    empty = []
    rows = []
    dens = Fraction(1)
    degs = degenerate_ts(m, n)
    for p in P:
        if p == 2: continue
        sv = survivors_mod_p(m, n, p)
        # какие из выживших — редукции известных вырожденных t
        degred = set()
        for d in degs:
            try:
                r = (d.numerator * pow(d.denominator, -1, p)) % p
                degred.add(r)
            except ValueError:
                degred.add('inf')
        extra = [x for x in sv if x not in degred]
        rows.append((p, len(sv), len(extra)))
        dens *= Fraction(len(sv), p+1)
        if len(sv) == 0: empty.append(p)
    for p, c, e in rows:
        print(f"  p={p:3d}: выжило {c:3d} из {p+1:3d} точек P^1(F_p); "
              f"не редукции вырожденных t: {e}")
    print()
    print(f"  ПРОСТЫХ С ПУСТЫМ ОБРАЗОМ (локальное препятствие): {empty if empty else 'НЕТ НИ ОДНОГО'}")
    print(f"  Произведение плотностей по всем нечётным p<200: {float(dens):.6e}")
    print()
    print("="*78)
    print("2. (15,8): 2-адический тест (a,b mod 2^N, точная арифметика Fraction)")
    print("="*78)
    for N in (3, 4, 5, 6):
        g = two_adic_survivors(m, n, N)
        print(f"  N={N}: пар (a,b) mod 2^{N} с пятью квадратами в Q_2: {len(g)} из {(1<<N)*(1<<N)}")
        if g and N == 4:
            print(f"        примеры: {g[:8]}")
    print()
    print("="*78)
    print("3. КОНТРОЛЬ: воспроизведение исключений Codex (13,8)@17 и (16,5)@41")
    print("="*78)
    for (mm, nn) in [(13,8), (16,5)]:
        firstempty = None
        for p in P:
            if p == 2: continue
            if len(survivors_mod_p(mm, nn, p)) == 0:
                firstempty = p; break
        print(f"  ({mm},{nn}): первое простое с пустым образом = {firstempty}")
    print()
    print("="*78)
    print("4. ВСЕ 127 ПАР m,n<=20, gcd=1: у кого есть локальное препятствие при p<200")
    print("="*78)
    pairs = [(a,b) for a in range(2,21) for b in range(1,a) if math.gcd(a,b)==1 and (a-b)%2==1]
    # G1 требует m>n>=1, gcd=1; s=(m^2+n^2)/2 полуцелое => m,n разной чётности
    print(f"  пар в списке: {len(pairs)}")
    blocked, survive = [], []
    for (mm,nn) in pairs:
        fe = None
        for p in P:
            if p == 2: continue
            if len(survivors_mod_p(mm,nn,p)) == 0:
                fe = p; break
        if fe: blocked.append(((mm,nn),fe))
        else:  survive.append((mm,nn))
    print(f"  С ЛОКАЛЬНЫМ ПРЕПЯТСТВИЕМ (p<200): {len(blocked)}")
    for pr, p in blocked: print(f"     {pr} -> p={p}")
    print(f"  БЕЗ ЛОКАЛЬНОГО ПРЕПЯТСТВИЯ (p<200): {len(survive)}")
    print("   ", survive)
