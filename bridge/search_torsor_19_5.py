#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
САМАЯ ОСТРАЯ ПОПЫТКА СЛОМАТЬ исключение (19,5).

Заявление Codex эквивалентно тому, что класс delta=(1,sigma,sigma) лежит в Ш(E/Q)[2]:
2-накрытие H_delta всюду локально разрешимо (проверено отдельно), но рациональных
точек не имеет. ОДНА рациональная точка на H_delta дала бы rank E(Q)>=2 и обрушила
бы весь аргумент.

Общий вид H_delta для delta=(1,sigma,sigma), sigma = бесквадратная часть s:
    w1^2 - sigma*w2^2 = (e2-e1) w0^2
    sigma*(w2^2 - w3^2) = (e3-e2) w0^2
После очистки знаменателей (w0 = q*v):
    w2^2 - w3^2 = K v^2            (решается ТОЧНО факторизацией K v^2 = u*w)
    w1^2 = sigma*w2^2 + C v^2
Перебираем v и ВСЕ разложения K v^2 = u*w — это не бокс-поиск, а точное решение
первого уравнения, поэтому глубина по (w2,w3) фактически неограничена.

ПОЗИТИВНЫЙ КОНТРОЛЬ: тот же код на (15,8), где точка ИЗВЕСТНА (t=1, u=17).
"""
from math import isqrt, gcd
import sys


def issq(x):
    if x < 0:
        return False
    r = isqrt(x)
    return r * r == x


def factorize(mm):
    f = {}
    x = mm
    p = 2
    while p * p <= x:
        while x % p == 0:
            f[p] = f.get(p, 0) + 1
            x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        f[x] = f.get(x, 0) + 1
    return f


def divisors_from(fd):
    ds = [1]
    for p, e in fd.items():
        nd = []
        pe = 1
        for _ in range(e + 1):
            nd.extend(x * pe for x in ds)
            pe *= p
        ds = nd
    return ds


def search(K, C, sigma, VMAX, label, verbose_hits=True):
    """ищет (v,w1,w2,w3): w2^2-w3^2 = K v^2 и w1^2 = sigma*w2^2 + C v^2"""
    KF = factorize(K)
    hits = []
    cnt = 0
    for v in range(1, VMAX + 1):
        tot = dict(KF)
        for p, e in factorize(v).items():
            tot[p] = tot.get(p, 0) + 2 * e
        Nn = K * v * v
        for u in divisors_from(tot):
            w = Nn // u
            if u > w or (u + w) % 2:
                continue
            w2 = (u + w) // 2
            w3 = (w - u) // 2
            cnt += 1
            val = sigma * w2 * w2 + C * v * v
            if val < 0 or not issq(val):
                continue
            w1 = isqrt(val)
            if gcd(gcd(w1, w2), gcd(w3, v)) != 1:
                pass
            hits.append((v, w1, w2, w3))
            if verbose_hits:
                print("   ТОЧКА %s: v=%d w1=%d w2=%d w3=%d" % (label, v, w1, w2, w3))
            if len(hits) > 4:
                return hits, cnt
    return hits, cnt


print("=== ПОЗИТИВНЫЙ КОНТРОЛЬ: (m,n)=(15,8), delta=(1,2,2), известна точка t=1,u=17 ===")
# s=289/2, sigma=2, e1=-2080800, e2=-14630625/2, e3=-591872
# w2^2-w3^2 = 13446881 v^2 ;  w1^2 = 2 w2^2 - 20938050 v^2   (w0 = 2v)
ctrl, ccnt = search(13446881, -20938050, 2, 3, "(15,8)")
print("   проверено пар (w2,w3): %d ; найдено точек: %d" % (ccnt, len(ctrl)))
ok = any((v, w1, w2, w3) == (1, 4080, 4335, 2312) for (v, w1, w2, w3) in ctrl)
print("   ожидаемая точка (v,w1,w2,w3)=(1,4080,4335,2312) найдена? %s" % ok)
print("   => СХЕМА ПОИСКА РАБОТАЕТ" if ok else "   => КОНТРОЛЬ ПРОВАЛЕН, поиску верить нельзя")
print()

print("=== ЦЕЛЬ: (m,n)=(19,5), delta=(1,193,193) ===")
# s=193 целое, sigma=193, w0=v.
# w2^2-w3^2 = (e3-e2)/193 * v^2 = 129696 v^2
# w1^2 = 193 w2^2 + (e2-e1) v^2 = 193 w2^2 - 23410128 v^2
VMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
hits, cnt = search(129696, -23410128, 193, VMAX, "(19,5)")
print("v до %d; проверено ТОЧНЫХ решений первого уравнения (w2,w3): %d" % (VMAX, cnt))
print("найдено рациональных точек H_(1,193,193): %d" % len(hits))
if hits:
    print("!!! ЗАЯВЛЕНИЕ CODEX ОПРОВЕРГНУТО:", hits[:3])
else:
    print("точек не найдено (это НЕ доказательство отсутствия, а согласованность)")
