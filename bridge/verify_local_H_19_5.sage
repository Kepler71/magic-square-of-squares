# verify_local_H_19_5.sage -- независимая перепроверка local_H_19_5.sage
load("/home/kep/magicKube/bridge/local_H_19_5.sage")

print()
print("#"*72)
print("НЕЗАВИСИМАЯ ВЕРИФИКАЦИЯ")
print("#"*72)

# ---------- A. проверка теста квадратичности в Q_p против Sage Qp ----------
print("\nA. is_square_in_Qp против Sage Qp(...).is_square()")
bad_cases = 0
set_random_seed(1)
for p in [2,3,5,7,11,13,19,193]:
    K = Qp(p, prec=60)
    for _ in range(400):
        x = QQ(ZZ.random_element(-10^6, 10^6))/ZZ.random_element(1, 1000)
        if x == 0: continue
        mine = is_square_in_Qp(x, p)
        theirs = K(x).is_square()
        if mine != theirs:
            bad_cases += 1
            print("   РАСХОЖДЕНИЕ p=%s x=%s mine=%s sage=%s" % (p,x,mine,theirs))
print("   расхождений:", bad_cases)

# ---------- B. проверка рекурсии на заведомо известных случаях ----------
print("\nB. exists_t_with_square на контрольных примерах")
PZ.<S> = ZZ[]
tests = [
    (PZ(-1), 3, False, "f=-1, p=3: -1 не квадрат в Q_3"),
    (PZ(-1), 5, True,  "f=-1, p=5: -1 квадрат в Q_5"),
    (PZ(3),  3, False, "f=3, p=3: v=1 нечётна"),
    (PZ(2),  2, False, "f=2, p=2: v=1 нечётна"),
    (PZ(17), 2, True,  "f=17, p=2: 17=1 mod 8"),
    (PZ(5),  2, False, "f=5, p=2: 5!=1 mod 8"),
    (PZ(S),  7, True,  "f=t, p=7: t=1"),
    (3*PZ(S)^2+3, 3, None, "f=3(t^2+1), p=3: см. ниже"),
]
for F, p, expect, desc in tests:
    r = exists_t_with_square(F, p, maxdepth=40)
    ok = "" if expect is None else ("OK" if r == expect else "!!! ОШИБКА")
    print("   %-45s -> %s  %s" % (desc, r, ok))
# 3(t^2+1) над Q_3: t^2+1 = 1,2,2 mod 3 -> никогда 0 mod 3, v_3=1 нечётна => False
print("   (ожидание для 3(t^2+1) над Q_3: False, т.к. t^2+1 не делится на 3)")

# ---------- C. рандомизированная сверка рекурсии с грубым перебором ----------
print("\nC. рандомизированная сверка: рекурсия vs грубый перебор t in {0..p^k-1}")
set_random_seed(7)
mismatch = 0
for trial in range(300):
    p = choice([2,3,5,7,11,13])
    deg = choice([2,3,4,5,6])
    F = PZ([ZZ.random_element(-40,40) for _ in range(deg+1)])
    if F.is_zero(): continue
    r = exists_t_with_square(F, p, maxdepth=25)
    k = 10 if p in (2,3) else 6
    brute = any(is_square_in_Qp(F(t0), p) for t0 in range(p^k))
    if r is True and not brute:
        mismatch += 1; print("   !!! рекурсия True, перебор не нашёл:", F, p)
    if r is False and brute:
        mismatch += 1; print("   !!! рекурсия False, перебор НАШЁЛ:", F, p)
print("   расхождений:", mismatch)

# ---------- D. ручная проверка сертификатов для (19,5) ----------
print("\nD. ручная проверка сертификатов для H(19,5)")
certs = {2: 0, 3: 0, 5: 1, 7: 0, 11: 1, 13: 1, 19: 1, 193: 112}
for p, t0 in certs.items():
    x = f(t0)
    e = ZZ(x).valuation(p)
    u = ZZ(x)/p^e
    if p == 2:
        cond = (u % 8 == 1)
        det = "u mod 8 = %s" % (u % 8)
    else:
        cond = (kronecker(ZZ(u), p) == 1)
        det = "kronecker(u,%s) = %s" % (p, kronecker(ZZ(u), p))
    print("   p=%-4s t=%-6s f(t)=%s" % (p, t0, x))
    print("        v_%s = %s (чётно: %s), %s  =>  квадрат в Q_%s: %s"
          % (p, e, e % 2 == 0, det, p, (e % 2 == 0) and cond))
    assert (e % 2 == 0) and cond, "сертификат p=%s НЕВЕРЕН" % p

# ---------- E. проверка через Qp: точка реально поднимается ----------
print("\nE. подъём Y в Q_p (проверка sqrt в Sage)")
for p, t0 in certs.items():
    K = Qp(p, prec=40)
    Y = K(f(t0)).sqrt()
    resid = K(Y^2 - f(t0))
    print("   p=%-4s t=%-6s Y = %s ...,  Y^2 - f(t) = %s" % (p, t0, str(Y)[:40], resid))
    assert resid == 0 or resid.valuation() > 30

# ---------- F. хорошая редукция: аккуратная проверка ----------
print("\nF. хорошая редукция и граница Вейля")
print("   bad =", bad)
for p in prime_range(3, 60):
    if p in bad: continue
    assert disc_f % p != 0 and lc % p != 0
print("   для всех нечётных p не из bad: p не делит ни disc(f), ни lc(f)  => модель")
print("   Y^2=f(t) гладка над Z_p, редукция -- гладкая кривая рода 2 над F_p.")
for p in [17, 19, 23, 29]:
    print("   p=%-3s : p+1-4*sqrt(p) = %.3f" % (p, p + 1 - 4*sqrt(float(p))))
print("   => при p>=17 (хорошая редукция) #H(F_p) > 0, гладкая точка поднимается по Гензелю.")

# ---------- G. точное число точек mod p для p < 17 (для полноты) ----------
print("\nG. #H(F_p) для малых p (аффинные + бесконечность, только гладкая модель при хорошей редукции)")
for p in [3,5,7,11,13]:
    Fp = GF(p); cnt = 0
    for t0 in range(p):
        v = Fp(f(t0))
        cnt += 2 if (v != 0 and v.is_square()) else (1 if v == 0 else 0)
    if Fp(lc) != 0 and Fp(lc).is_square(): cnt += 2
    elif Fp(lc) == 0: cnt += 1
    print("   p=%-3s : #точек на модели = %s" % (p, cnt))

print("\nВЕРИФИКАЦИЯ ЗАВЕРШЕНА")
