# -*- coding: utf-8 -*-
# Независимая ВТОРАЯ реализация проверки локального препятствия (16,5) mod 41.
# Ничего не берёт из питоновского скрипта: работает прямо в GF(p),
# перебирает ВСЕ пары (a,b) != (0,0) в F_p^2, а не проективных представителей.
# Запуск: sage /home/kep/magicKube/bridge/local_16_5_p41_скептик_sage.sage

def five_cells_GF(m, n, a, b, F):
    two = F(2)
    s = (F(m)^2 + F(n)^2) / two
    F0 = F(m)^2*b^2 + F(n)^2*a^2
    F4 = s*(a^2 + b^2)
    F8 = F(n)^2*b^2 + F(m)^2*a^2
    L  = F4 - two*F(m)*F(n)*a*b
    U  = F4 + two*F(m)*F(n)*a*b
    return [F0, F4, F8, L, U]

def scan_all_affine_pairs(m, n, p, zero_ok=True):
    """перебор ВСЕХ (a,b) в F_p^2 \\ {(0,0)} — без понятия проективного представителя"""
    F = GF(p)
    good = []
    for a in F:
        for b in F:
            if a == 0 and b == 0:
                continue
            vals = five_cells_GF(m, n, a, b, F)
            ok = all(v.is_square() and (zero_ok or v != 0) for v in vals)
            if ok:
                good.append((a, b, vals))
    return good

def scan_projective(m, n, p, idx=(0,1,2,3,4), zero_ok=True):
    F = GF(p)
    P = ProjectiveSpace(1, F)
    good = []
    for pt in P.rational_points():
        a, b = pt[0], pt[1]
        vals = five_cells_GF(m, n, a, b, F)
        if all(vals[i].is_square() and (zero_ok or vals[i] != 0) for i in idx):
            good.append((a, b, vals))
    return good

print("="*100)
print("S1. Вторая независимая реализация: (16,5), p=41, перебор всех 1680 пар (a,b) в F_41^2")
print("="*100)
p = 41
good = scan_all_affine_pairs(16, 5, p)
print("  пар (a,b) != (0,0), где все пять клеток — квадраты в F_41:", len(good))
print("  тот же счёт через ProjectiveSpace(1,GF(41)).rational_points():",
      len(scan_projective(16, 5, p)))
print("  число точек P^1(F_41):", len(ProjectiveSpace(1, GF(p)).rational_points()))

print()
print("  прошедшие только три клетки C7 (проективно):")
for a, b, v in scan_projective(16, 5, p, idx=(0,1,2)):
    print(f"     (a:b)=({a}:{b})  (F0,F4,F8)={tuple(v[:3])}  (L,U)={tuple(v[3:])}  "
          f"L кв? {v[3].is_square()}  U кв? {v[4].is_square()}")

print()
print("  тот же счёт для перестановки (m,n)=(5,16):", len(scan_projective(5, 16, p)))
print("  тот же счёт для масштабирования (m,n)=(32,10) [то же отношение]:",
      len(scan_projective(32, 10, p)))
print("  тот же счёт для (m,n)=(16+41, 5) [та же редукция]:",
      len(scan_projective(57, 5, p)))

print()
print("="*100)
print("S2. Все НЕЧЁТНЫЕ простые p < 300, включая делители 2mn(m^2-n^2)(m^2+n^2)")
print("    (аргумент редукции не требует хорошей редукции — только p нечётно)")
print("="*100)
bad_div = 2*16*5*(16^2-5^2)*(16^2+5^2)
killers = []
for q in prime_range(3, 300):
    sv = scan_projective(16, 5, q)
    tag = " [делит 2mn(m^2-n^2)(m^2+n^2)]" if bad_div % q == 0 else ""
    if len(sv) == 0:
        killers.append(q)
        print(f"   p={q:4d}: ПУСТО — препятствие{tag}")
print("   все простые с пустым P^1(F_p):", killers)
print("   наименьшее:", killers[0] if killers else "нет")

print()
print("="*100)
print("S3. Контроль: кривая C7 (три клетки) для (16,5) — есть ли препятствие на семи клетках?")
print("="*100)
k3 = [q for q in prime_range(3, 300) if len(scan_projective(16, 5, q, idx=(0,1,2))) == 0]
print("   простые с пустым C7:", k3 or "нет — C7 локально всюду непусто, препятствие даёт ТОЛЬКО пара L,U")

print()
print("="*100)
print("S4. Контроль на парах с известными рациональными точками (тест не должен их убивать)")
print("="*100)
for (m, n) in [(7,1), (15,8), (4,3), (13,8)]:
    ks = [q for q in prime_range(3, 300) if len(scan_projective(m, n, q)) == 0]
    print(f"   (m,n)=({m},{n}): простые с пустым P^1(F_p): {ks or 'нет'}")

print()
print("="*100)
print("S5. Проверка того, что у (16,5) нет рациональной точки малой высоты на пяти клетках")
print("    (глобальный контроль; отсутствие находки НЕ доказательство)")
print("="*100)
found = []
H = 120
for b in range(1, H+1):
    for a in range(-H, H+1):
        if gcd(abs(a), b) != 1:
            continue
        m, n = 16, 5
        F0 = m^2*b^2 + n^2*a^2
        F4 = QQ(m^2+n^2)*(a^2+b^2)/2
        F8 = n^2*b^2 + m^2*a^2
        L  = F4 - 2*m*n*a*b
        U  = F4 + 2*m*n*a*b
        if all(QQ(x) >= 0 and QQ(x).is_square() for x in [F0, F4, F8, L, U]):
            found.append((a, b))
print(f"   t=a/b, |a|,b <= {H}: найдено {len(found)} -> {found}")

print()
print("="*100)
print("S6. Что именно ломается: распределение (L,U) по четырём выжившим C7-точкам")
print("="*100)
F = GF(41)
for t in [4, 10, 31, 37]:
    v = five_cells_GF(16, 5, F(t), F(1), F)
    print(f"   t={t:2d}: F0={v[0]} (кв {v[0].is_square()}), F4={v[1]} (кв {v[1].is_square()}), "
          f"F8={v[2]} (кв {v[2].is_square()}), L={v[3]} (кв {v[3].is_square()}), "
          f"U={v[4]} (кв {v[4].is_square()}), L*U={v[3]*v[4]} (кв {(v[3]*v[4]).is_square()})")
v = five_cells_GF(16, 5, F(1), F(0), F)
print(f"   inf : F0={v[0]} (кв {v[0].is_square()}), F4={v[1]} (кв {v[1].is_square()}), "
      f"F8={v[2]} (кв {v[2].is_square()})")
