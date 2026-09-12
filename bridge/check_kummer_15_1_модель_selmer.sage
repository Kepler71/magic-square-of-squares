# -*- coding: utf-8 -*-
# Спутник к check_kummer_15_1_модель.sage
# ПЕРЕКРЁСТНАЯ АТАКА НА СОГЛАСОВАННОСТЬ:
#   Если C всюду локально разрешима (заявление Codex §6), то требуемый класс
#   (1,113,113) ОБЯЗАН лежать в 2-Селмеровой группе E.  Тогда
#       dim Sel_2 = r + 2 + dim Sha[2] >= 3 + 1,  а Sha[2] чётномерна => >= 5.
#   Если же Sage даст dim Sel_2 = 3, то Sha[2]=0, класс НЕ в Sel_2, и тогда
#   заявление о всюду локальной разрешимости C ЛОЖНО (исключение при этом
#   осталось бы верным, но «глобальность барьера» — нет).
# Плюс: собственная проверка локальной разрешимости C.

import sys

def hdr(s):
    print("\n" + "=" * 78); print(s); print("=" * 78); sys.stdout.flush()

FAIL = []
def check(name, cond, extra=""):
    ok = bool(cond)
    print(("  [OK]   " if ok else "  [ПРОВАЛ] ") + name + (("  " + extra) if extra else ""))
    sys.stdout.flush()
    if not ok: FAIL.append(name)
    return ok

m0, n0 = ZZ(15), ZZ(1)
s0 = QQ(m0^2 + n0^2)/2
b0 = s0*m0^2*n0^2
E1 = -b0; E2 = -s0*m0^4; E3 = -s0*n0^4
Ecurve = EllipticCurve([0, -(E1+E2+E3), 0, E1*E2+E1*E3+E2*E3, -E1*E2*E3])

hdr("1. 2-Селмерова группа и согласованность с локальной разрешимостью C")
print("  E =", Ecurve)
sys.stdout.flush()
sel = None
try:
    sel = Ecurve.selmer_rank()
    print("  Sage selmer_rank() =", sel)
except Exception as ex:
    print("  selmer_rank:", type(ex).__name__, str(ex)[:200])
try:
    print("  rank_bound() =", Ecurve.rank_bound())
except Exception as ex:
    print("  rank_bound:", ex)
try:
    pr = pari(Ecurve.minimal_model().ainvs()).ellinit().ellrank()
    print("  PARI ellrank =", pr)
except Exception as ex:
    print("  PARI:", ex)

if sel is not None:
    print("\n  Логика: dim Sel_2 = rank + 2 + dim Sha[2]  (полное 2-кручение).")
    print(f"  при rank = 1: dim Sha[2] = {sel} - 3 = {sel - 3}")
    if sel == 3:
        print("  !!! Sha[2] = 0 => класс (1,113,113) НЕ в Sel_2 =>")
        print("  !!! C НЕ всюду локально разрешима, §6 отчёта Codex ложен.")
    else:
        check("dim Sha[2] > 0 (совместимо с всюду локальной разрешимостью C)", sel - 3 > 0)
        check("dim Sha[2] чётна (пары Касселса)", (sel - 3) % 2 == 0, f"dim Sha[2] = {sel-3}")

# ----------------------------------------------------------------------------
hdr("2. Собственная проверка: локальная разрешимость C (ищу t сам)")
F0 = lambda t: m0^2 + n0^2*t^2
F4 = lambda t: s0*(1 + t^2)
F8 = lambda t: n0^2 + m0^2*t^2

def is_sq_qp(a, p):
    """a in Q*, проверка: a — ненулевой квадрат в Q_p."""
    a = QQ(a)
    if a == 0: return False
    v = a.valuation(p)
    if v % 2 != 0: return False
    u = a / QQ(p)^v
    num, den = u.numerator(), u.denominator()
    if p == 2:
        # u — 2-адическая единица; квадрат <=> u = 1 mod 8
        return (num * inverse_mod(den % 8, 8)) % 8 == 1
    return kronecker(num * den, p) == 1

badp = [2, 3, 5, 7, 113]
allp = badp + [p for p in primes(200) if p not in badp]
found = {}
missing = []
for p in allp:
    hit = None
    # ищем t = a/c с малыми a,c и достаточной p-адической «глубиной»
    cands = [QQ(a) for a in range(-40, 41)]
    cands += [QQ(a)/QQ(c) for a in range(-20, 21) for c in range(1, 21) if gcd(a, c) == 1]
    cands += [QQ(p)^k * QQ(a) for k in range(-3, 4) for a in range(-9, 10) if a != 0]
    for t in cands:
        f0, f4, f8 = F0(t), F4(t), F8(t)
        if f0 != 0 and f4 != 0 and f8 != 0 and \
           is_sq_qp(f0, p) and is_sq_qp(f4, p) and is_sq_qp(f8, p):
            hit = t; break
    if hit is not None:
        found[p] = hit
    else:
        missing.append(p)
print("  найдены явные t для p:", sorted(found.keys())[:40], "..." if len(found) > 40 else "")
print("  примеры:", {p: found[p] for p in badp if p in found})
if missing:
    print("  НЕ найдено t поиском для p =", missing)
check("для всех плохих p={2,3,5,7,113} найден явный t", all(p in found for p in badp))
# R
check("C(R) != пусто (t=0 даёт положительные F)", F0(0) > 0 and F4(0) > 0 and F8(0) > 0)

# для больших хороших p — счёт точек на C над F_p напрямую
hdr("3. Хорошие p: прямой подсчёт точек C(F_p) (не Хассе-Вейль на слово)")
def count_C_Fp(p):
    if p in badp: return None
    Fp = GF(p)
    m_, n_, s_ = Fp(m0), Fp(n0), Fp(s0)
    cnt = 0
    for t in Fp:
        v0 = m_^2 + n_^2*t^2; v4 = s_*(1+t^2); v8 = n_^2 + m_^2*t^2
        if v0 == 0 or v4 == 0 or v8 == 0: continue
        if v0.is_square() and v4.is_square() and v8.is_square():
            cnt += 8
    return cnt
worst = []
for p in primes(101, 400):
    if p in badp: continue
    c = count_C_Fp(p)
    if c == 0: worst.append(p)
check("C(F_p) непусто для всех хороших 101<=p<400 (гладкие точки поднимаются по Гензелю)",
      not worst, f"пустые: {worst}")
print("  для p>=101 также годится Хассе-Вейль: p+1-10*sqrt(p) > 0 при p >= 101")
check("Хассе-Вейль даёт >0 при p>=101", 101 + 1 - 10*sqrt(101.0) > 0,
      f"101+1-10*sqrt(101) = {101+1-10*sqrt(101.0):.3f}")

# ----------------------------------------------------------------------------
hdr("4. Торсор D для класса (1,113,113): проверяю его вывод сам")
A = E2 - E1; B = E3 - E1
print(f"  A = e2-e1 = {A},  B = e3-e1 = {B}")
print(f"  D: R1^2 - 113 R2^2 = {A} W^2 ;  R1^2 - 113 R3^2 = {B} W^2")
check("A = -5695200 (совпадает с Codex)", A == -5695200)
check("B = 25312  (совпадает с Codex)", B == 25312)
# отображение C -> D: (R1,R2,R3,W) = (m*u4, m*u0, u8, 1)? проверим символьно
Rsym = PolynomialRing(QQ, ['t']); tv = Rsym.gen()
f0 = m0^2 + n0^2*tv^2; f4 = s0*(1+tv^2); f8 = n0^2 + m0^2*tv^2
# кандидат: R1^2 = X - e1 = (m n)^2 u4^2 = 225 f4 ; R2^2 = (X-e2)/113 = 225 f0 ; R3^2 = (X-e3)/113 = f8
lhs1 = Rsym(225*f4) - 113*Rsym(225*f0)
lhs2 = Rsym(225*f4) - 113*Rsym(f8)
check("R1=15u4, R2=15u0, R3=u8, W=1 удовлетворяет 1-му уравнению D", lhs1 == Rsym(A))
check("R1=15u4, R2=15u0, R3=u8, W=1 удовлетворяет 2-му уравнению D", lhs2 == Rsym(B))
print("  => C -> D задано (t,u0,u4,u8) |-> (15u4, 15u0, u8, 1). Подтверждено символьно.")
# W=0 ветвь
check("при W=0: R1^2=113R2^2 не имеет нетривиальных рациональных решений (113 не квадрат)",
      not QQ(113).is_square())

hdr("ИТОГ спутника")
if FAIL:
    print("ПРОВАЛЕННЫЕ ПРОВЕРКИ:")
    for f in FAIL: print("   -", f)
else:
    print("Все проверки спутника пройдены.")
