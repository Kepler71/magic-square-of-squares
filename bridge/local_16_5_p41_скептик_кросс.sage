# -*- coding: utf-8 -*-
# СТАТУС: независимая перепроверка (другая реализация, другой аппарат)
# ДЛЯ: сверка с local_16_5_p41_скептик.py
# ИТОГ: точки системы девяти клеток над F_41 для (16,5) перечисляются ПЕРЕБОРОМ РЕШЕНИЙ, а не тестом QR
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЕНО: та же процедура на контроле (15,8), где решение над Q известно
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: сомневаетесь в моей же реализации символа Лежандра
#
# Запуск: sage /home/kep/magicKube/bridge/local_16_5_p41_скептик_кросс.sage

def brute_points(m, n, p):
    """Перебор ВСЕХ (a:b) in P^1(F_p) и всех (u0,u4,u8,uL,uU) in F_p^5 —
       ищем настоящие решения системы, не пользуясь никакими тестами на квадратичность."""
    F = GF(p)
    sq = {}                      # значение -> есть ли корень
    for y in F:
        sq[y*y] = True
    m2, n2 = F(m)^2, F(n)^2
    s = (F(m)^2 + F(n)^2) / F(2)
    mn2 = 2*F(m)*F(n)
    pts = [(F(t), F(1)) for t in range(p)] + [(F(1), F(0))]
    good3, good5 = [], []
    for (a, b) in pts:
        F0 = m2*b^2 + n2*a^2
        F4 = s*(a^2 + b^2)
        F8 = n2*b^2 + m2*a^2
        L  = F4 - mn2*a*b
        U  = F4 + mn2*a*b
        # честный перебор корней, без символа Лежандра
        r3 = all(any(y*y == v for y in F) for v in (F0, F4, F8))
        r5 = r3 and all(any(y*y == v for y in F) for v in (L, U))
        lab = str(a) if b == 1 else 'oo'
        if r3: good3.append(lab)
        if r5: good5.append(lab)
    return good3, good5, len(pts)

def nine_cells_mod_p(m, n, p, a, b):
    """все ДЕВЯТЬ однородных клеток mod p — проверяем, что четыре «автоматические» действительно квадраты"""
    F = GF(p)
    a, b = F(a), F(b)
    s = (F(m)^2 + F(n)^2)/F(2)
    c4 = s*(a^2+b^2)
    cells = [F(m)^2*b^2 + F(n)^2*a^2, (F(m)*a+F(n)*b)^2, c4 - 2*F(m)*F(n)*a*b,
             (F(m)*a-F(n)*b)^2, c4, (F(m)*b+F(n)*a)^2,
             c4 + 2*F(m)*F(n)*a*b, (F(m)*b-F(n)*a)^2, F(n)^2*b^2 + F(m)^2*a^2]
    return cells

print("="*90)
print("Кросс-проверка перебором РЕШЕНИЙ над F_p (без символа Лежандра)")
print("="*90)
for (m, n, p) in [(16,5,41), (5,16,41), (13,8,17), (15,8,41), (15,8,17)]:
    g3, g5, N = brute_points(m, n, p)
    print(f"  (m,n)=({m},{n}), p={p}: проверено {N} проективных точек; "
          f"три клетки -> {g3}; ПЯТЬ клеток -> {g5 if g5 else 'НЕТ НИ ОДНОЙ'}")

print()
print("="*90)
print("Контроль магичности и автоматических квадратов mod 41 в каждой точке")
print("="*90)
m, n, p = 16, 5, 41
F = GF(p)
bad_magic, bad_auto = [], []
for (a, b) in [(t,1) for t in range(p)] + [(1,0)]:
    c = nine_cells_mod_p(m, n, p, a, b)
    tot = 3*c[4]
    lines = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
    if any(c[i]+c[j]+c[k] != tot for (i,j,k) in lines):
        bad_magic.append((a,b))
    if any(not any(y*y == c[i] for y in F) for i in (1,3,5,7)):
        bad_auto.append((a,b))
print(f"  точек с нарушенной магичностью mod 41: {len(bad_magic)} (должно быть 0)")
print(f"  точек, где какая-то из клеток c1,c3,c5,c7 НЕ квадрат mod 41: {len(bad_auto)} (должно быть 0)")

print()
print("="*90)
print("Полная система как схема: перебор всех решений (a,b,u0,u2,u4,u6,u8) над F_41")
print("="*90)
# прямой перебор по u-переменным: ищем (a:b) и корни для всех пяти клеток одновременно
cnt = 0
for (a, b) in [(t,1) for t in range(p)] + [(1,0)]:
    c = nine_cells_mod_p(m, n, p, a, b)
    free = [c[0], c[2], c[4], c[6], c[8]]
    roots = [[y for y in F if y*y == v] for v in free]
    cnt += prod([len(r) for r in roots])
print(f"  общее число решений (a:b, u0,u2,u4,u6,u8) над F_41 для (16,5): {cnt}")
g3, g5, _ = brute_points(15, 8, 41)
cnt2 = 0
for (a, b) in [(t,1) for t in range(p)] + [(1,0)]:
    c = nine_cells_mod_p(15, 8, p, a, b)
    free = [c[0], c[2], c[4], c[6], c[8]]
    roots = [[y for y in F if y*y == v] for v in free]
    cnt2 += prod([len(r) for r in roots])
print(f"  то же для КОНТРОЛЬНОЙ (15,8) над F_41: {cnt2} (обязано быть > 0)")
