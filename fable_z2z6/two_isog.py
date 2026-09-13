# -*- coding: utf-8 -*-
"""Точные группы Селмера 2-изогении для E(a,b) из каждой из трёх точек порядка 2.
Локальная разрешимость квартик — точным решателем Claude criterion_proof/selmer_exact.py (locally_solvable).
Свой драйвер (Fable), 14.09.2026."""
import sys, os
sys.path.insert(0, '/home/kep/magicKube/criterion_proof')
from selmer_exact import locally_solvable, primes_of
from itertools import combinations
from math import log2

def squarefree_divisors_signed(n, positive_only=False):
    P = primes_of(n)
    out = []
    for r in range(len(P)+1):
        for sub in combinations(P, r):
            d = 1
            for q in sub: d *= q
            out.append(d)
            if not positive_only: out.append(-d)
    return out

def real_solvable(c0, c2, c4):
    """w^2 = c4 z^4 + c2 z^2 + c0 имеет вещественную точку?"""
    if c0 >= 0 or c4 > 0: return True
    return c2 > 0 and c2*c2 - 4*c4*c0 >= 0

def isog_selmer(A, B, D):
    """E: y^2 = x(x^2+Ax+B), E': Y^2 = X(X^2-2AX+D), D=A^2-4B=(h-k)^2 — квадрат.
    Возвращает (Sphi, Sphihat)."""
    assert A*A - 4*B == D
    hk = int(round(D**0.5))
    while hk*hk < D: hk += 1
    while hk*hk > D: hk -= 1
    assert hk*hk == D
    bad = sorted(set([2] + primes_of(B) + primes_of(hk)))
    S = []
    for d in squarefree_divisors_signed(B):
        c = [B//d, 0, A, 0, d]
        if not real_solvable(c[0], c[2], c[4]): continue
        if all(locally_solvable(c, p) for p in bad): S.append(d)
    T = []
    for d in squarefree_divisors_signed(hk):
        c = [D//d, 0, -2*A, 0, d]
        if not real_solvable(c[0], c[2], c[4]): continue
        if all(locally_solvable(c, p) for p in bad): T.append(d)
    return S, T

def roots_data(a, b):
    M = b**3*(2*a+b); N = a**3*(a+2*b)
    e = [0, -M, -N]
    out = []
    for i in range(3):
        j, k = [t for t in range(3) if t != i]
        h = e[i]-e[j]; kk = e[i]-e[k]
        out.append((h+kk, h*kk, (h-kk)**2))
    return out

def two_isog_bounds(a, b):
    res = []
    for (A, B, D) in roots_data(a, b):
        S, T = isog_selmer(A, B, D)
        res.append(dict(A=A, B=B, S=S, T=T, bound=int(round(log2(len(S)*len(T))))-2))
    return res

if __name__ == "__main__":
    a, b = int(sys.argv[1]), int(sys.argv[2])
    for i, r in enumerate(two_isog_bounds(a, b)):
        print(i, r['bound'], r['S'], r['T'])
