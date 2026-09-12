# -*- coding: utf-8 -*-
# СТАТУС: независимая перепроверка §4 скептика-III другим движком (Sage GF(p), is_square)
# ДЛЯ: пользователь, Codex
# ИТОГ: см. вывод
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЯЕТ: local_16_5_p41_скептик_III.py, §4 и §6
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужна вторая реализация того же перебора
#
# Запуск: sage /home/kep/magicKube/bridge/local_16_5_p41_скептик_III_sage.sage

def survivors(m, n, p, allow_zero=True, which=(0,1,2,3,4)):
    F = GF(p)
    s = F(m^2 + n^2) / F(2)
    out = []
    # все точки P^1(F_p) через саму Sage, а не через свой список
    for Pt in ProjectiveSpace(1, F).rational_points():
        a, b = F(Pt[0]), F(Pt[1])
        cells = [F(m)^2*b^2 + F(n)^2*a^2,
                 s*(a^2 + b^2),
                 F(n)^2*b^2 + F(m)^2*a^2,
                 s*(a^2 + b^2) - 2*F(m)*F(n)*a*b,
                 s*(a^2 + b^2) + 2*F(m)*F(n)*a*b]
        ok = True
        for i in which:
            c = cells[i]
            if c == 0:
                if not allow_zero:
                    ok = False; break
            elif not c.is_square():
                ok = False; break
        if ok:
            out.append((Pt, [ZZ(c) for c in cells]))
    return out

print("="*96)
print("Вторая реализация: Sage, ProjectiveSpace(1,GF(p)).rational_points(), c.is_square()")
print("="*96)

for (m, n, p) in [(16,5,41), (5,16,41), (13,8,17), (15,8,41), (15,8,17), (7,1,41), (11,4,41)]:
    sv5 = survivors(m, n, p)
    sv3 = survivors(m, n, p, which=(0,1,2))
    print(f"  (m,n)=({m},{n}), p={p}:  выживших на 5 клетках = {len(sv5)};  на 3 клетках (C7) = {len(sv3)}")
    if sv5:
        for Pt, c in sv5:
            print(f"       ВЫЖИЛА {Pt} -> {c}")
    if len(sv3) and len(sv3) <= 6:
        print(f"       (C7-выжившие: {[str(Pt) for Pt,_ in sv3]})")

print()
print("Число точек P^1(F_41), перечисленных Sage:", len(ProjectiveSpace(1, GF(41)).rational_points()))
print("Число точек P^1(F_17), перечисленных Sage:", len(ProjectiveSpace(1, GF(17)).rational_points()))

print()
print("Символическая сверка клеток (Sage, дробно-рациональные функции):")
R.<M,N,T> = QQ[]
Fr = Frac(R)
C = (M^2+N^2)*(1+T^2)/2
x = ((M-N*T)^2 - C + (M+N*T)^2 - C)/2
y = ((M+N*T)^2 - (M-N*T)^2)/2
cells = {'F0': C+x, 'F8': C-x, 'U': C+y, 'L': C-y, 'F4': C,
         'q1': C-x-y, 'q2': C-x+y, 'q3': C+x-y, 'q4': C+x+y}
print("   x =", Fr(x), "   y =", Fr(y))
for k in ['F0','F4','F8','L','U']:
    print(f"   {k} = {Fr(cells[k])}")
print("   четыре автоматических:")
for k in ['q1','q2','q3','q4']:
    print(f"   {k} = {factor(R(cells[k]))}")
codexF0 = M^2 + N^2*T^2
codexF8 = N^2 + M^2*T^2
codexL  = (M^2+N^2)*(1+T^2)/2 - 2*M*N*T
codexU  = (M^2+N^2)*(1+T^2)/2 + 2*M*N*T
print("   сверка с Codex: F0", Fr(cells['F0']-codexF0)==0,
      " F8", Fr(cells['F8']-codexF8)==0,
      " L", Fr(cells['L']-codexL)==0,
      " U", Fr(cells['U']-codexU)==0)
