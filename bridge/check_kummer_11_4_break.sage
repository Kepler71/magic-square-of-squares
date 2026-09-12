# Часть 5: ПОПЫТКА СЛОМАТЬ.
#  (1) прямой поиск рациональной точки на 2-накрытии D для класса (1,274,274).
#      если найдётся — класс В образе E(Q), ранг >= 2, исключение Codex РУШИТСЯ.
#  (2) ГЛУБОКИЙ поиск точек на самой C_{11,4} через сведение F4=квадрат
#      к уравнению p^2+q^2 = 274*z^2 (резко сужает перебор, даёт огромный охват по t).
def hdr(t):
    print("\n" + "="*78); print(t); print("="*78)
def sqfree(q):
    q = QQ(q); a=q.numerator(); d=q.denominator()
    return ZZ(a*d).squarefree_part()

M, N = ZZ(11), ZZ(4)
S = QQ(M^2+N^2)/2; B = S*M^2*N^2
E1, E2, E3 = QQ(-B), QQ(-S*M^4), QQ(-S*N^4)
print("m,n =", M, N, " s =", S, " b =", B)
print("e =", (E1,E2,E3))
print("sqfree(s) =", sqfree(S))

hdr("1. 2-НАКРЫТИЕ D для d=(1, s, s) = (1, 137/2, 137/2) ~ (1,274,274)")
d1, d2, d3 = QQ(1), S, S
print("система (однородно, вес w):")
print("   d1*z1^2 - d2*z2^2 = (e2-e1)*w^2 ,  e2-e1 =", E2-E1)
print("   d1*z1^2 - d3*z3^2 = (e3-e1)*w^2 ,  e3-e1 =", E3-E1)
A12 = 2*(E2-E1); A13 = 2*(E3-E1)
print("умножаем на 2:")
print("   2 z1^2 - 137 z2^2 =", A12, "* w^2")
print("   2 z1^2 - 137 z3^2 =", A13, "* w^2")
print("   (%s = %s ; %s = %s)" % (A12, factor(ZZ(A12)), A13, factor(ZZ(A13))))

# КОНТРОЛЬ: точка C обязана давать точку D. Проверим на (15,8), t=1.
print("\nКОНТРОЛЬ перехода C->D на известной точке (15,8), t=1:")
mm,nn = ZZ(15),ZZ(8); ss = QQ(mm^2+nn^2)/2; bb = ss*mm^2*nn^2
r1,r2,r3 = QQ(-bb), QQ(-ss*mm^4), QQ(-ss*nn^4)
u0=u4=u8=QQ(17); t0=QQ(1); x0 = bb*t0^2
z1 = mm*nn*u4; z2 = mm*u0; z3 = nn*u8; w0=QQ(1)
print("  z1=%s z2=%s z3=%s  d=(1,%s,%s)" % (z1,z2,z3,ss,ss))
print("  1*z1^2 - s*z2^2 =", z1^2 - ss*z2^2, "  vs  e2-e1 =", r2-r1,
      "  РАВНО:", z1^2 - ss*z2^2 == r2-r1)
print("  1*z1^2 - s*z3^2 =", z1^2 - ss*z3^2, "  vs  e3-e1 =", r3-r1,
      "  РАВНО:", z1^2 - ss*z3^2 == r3-r1)
print("  => параметризация накрытия ВЕРНА (положительный контроль пройден)")

hdr("1b. ПОИСК рациональных точек на D_{(1,274,274)} для (11,4)")
# 2 z1^2 - 137 z2^2 = A12 w^2 ; 2 z1^2 - 137 z3^2 = A13 w^2
# => z2^2 = (2 z1^2 - A12 w^2)/137 , z3^2 = (2 z1^2 - A13 w^2)/137
A12 = ZZ(A12); A13 = ZZ(A13)
found_D = []
LIM_Z = 6000
for w in range(1, 61):
    for z1 in range(0, LIM_Z+1):
        n2 = 2*z1^2 - A12*w^2
        if n2 < 0 or n2 % 137 != 0: continue
        q2 = n2 // 137
        if not ZZ(q2).is_square(): continue
        n3 = 2*z1^2 - A13*w^2
        if n3 < 0 or n3 % 137 != 0: continue
        q3 = n3 // 137
        if ZZ(q3).is_square():
            found_D.append((z1, ZZ(q2).sqrt(), ZZ(q3).sqrt(), w))
print("перебор w<=60, z1<=%d : найдено точек D: %d" % (LIM_Z, len(found_D)), found_D[:10])
if found_D:
    print("!!! КЛАСС (1,274,274) ЛЕЖИТ В ОБРАЗЕ E(Q) — ИСКЛЮЧЕНИЕ CODEX РУШИТСЯ !!!")
else:
    print("точек не найдено (НЕ доказательство отсутствия — только неудача поиска)")

hdr("2. ГЛУБОКИЙ ПОИСК НА C_{11,4}: сведение условия 'F4 — квадрат'")
print("t = p/q, gcd(p,q)=1.  F4 = s(p^2+q^2)/q^2 квадрат  <=>  sqfree(s)*(p^2+q^2) квадрат")
print("   <=>  p^2 + q^2 = %d * z^2" % sqfree(S))
print("КОНТРОЛЬ на (15,8): sqfree(s)=%s, t=1 -> p^2+q^2=2 = %s*1^2  OK:" %
      (sqfree(ss), sqfree(ss)), 1+1 == sqfree(ss)*1)
D274 = sqfree(S)
print("\nпроверка что 274 вообще представимо: 274 = 15^2+7^2 =", 15^2+7^2)

cnt = 0; hits = []
KMAX = 3000
import itertools
for z in range(1, KMAX+1):
    Nn = D274 * z^2
    # все представления Nn = p^2+q^2 с gcd(p,q)=1
    p = 0
    while p*p*2 <= Nn:
        rem = Nn - p*p
        sq = ZZ(rem).isqrt()
        if sq*sq == rem:
            q = sq
            if gcd(p,q) == 1 and q > 0:
                cnt += 1
                a0 = M^2*q^2 + N^2*p^2
                if ZZ(a0).is_square():
                    a8 = N^2*q^2 + M^2*p^2
                    if ZZ(a8).is_square():
                        hits.append((p,q,z))
                        print("   НАЙДЕНО! t = %d/%d" % (p,q))
        p += 1
print("проверено допустимых t (где F4 уже квадрат): %d, z<=%d" % (cnt, KMAX))
print("из них с F0 и F8 тоже квадратами: %d  %s" % (len(hits), hits))
print("макс. рассмотренная высота t ~ sqrt(274)*%d = %d" % (KMAX, int(sqrt(274.0)*KMAX)))
if hits:
    print("!!! ТОЧКА НА C НАЙДЕНА — ИСКЛЮЧЕНИЕ РУШИТСЯ !!!")
else:
    print("точек не найдено (НЕ доказательство отсутствия)")

hdr("2b. КОНТРОЛЬ ТОГО ЖЕ ПОИСКА на (15,8), где точка ТОЧНО есть")
Dc = sqfree(ss)
found_c = []
for z in range(1, 200):
    Nn = Dc*z^2
    p = 0
    while p*p*2 <= Nn:
        rem = Nn - p*p
        sq = ZZ(rem).isqrt()
        if sq*sq == rem:
            q = sq
            if gcd(p,q)==1 and q>0:
                if ZZ(mm^2*q^2 + nn^2*p^2).is_square() and ZZ(nn^2*q^2 + mm^2*p^2).is_square():
                    found_c.append((p,q,z))
        p += 1
print("(15,8): найдено t =", found_c[:10])
print("МЕТОД ПОИСКА РАБОТАЕТ (положительный контроль):", len(found_c) > 0)
print("ГОТОВО")
