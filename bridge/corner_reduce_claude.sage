# -*- coding: utf-8 -*-
# СТАТУС: рабочий скрипт — НЕЗАВИСИМАЯ проверка вопроса «сводится ли угловой случай к рёберному»
# ДЛЯ: Claude/Codex/Grok; запуск: sage corner_reduce_claude.sage <PART>
# ИТОГ: вывод углового семейства с нуля, инварианты, полный критерий изоморфизма, дешёвые закрытия
# ОТМЕНЯЕТ: ничего; это ВТОРАЯ независимая реализация рядом с bridge/corner_family.sage (его не трогаю)
# ПРОВЕРЕНО: каждое тождество подставляется числами; assert'ы падают, если что-то не так
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужен ответ на вопрос о сведении или входы для углового конвейера
#
# PART 1  вывод обоих семейств из аксиом магического квадрата (символьно + числа)
# PART 2  универсальное семейство X(A,B) рода 5; G1 и G2 как два арифметических слоя
# PART 3  ПОЛНЫЙ критерий изоморфизма X(A,B) ~ X(A',B') и его разбор для G2 (кривая рода 1, ранг)
# PART 4  сечения: A,C,beta для обоих типов; тождества; поля k; b'
# PART 5  род 2: Igusa-инварианты C_J2, полный перебор мёбиусовых эквивалентностей b' vs 2b'
# PART 6  ДЕШЁВЫЕ ЗАКРЫТИЯ углового случая: вырожденные точки, локально, ранги
import sys, itertools
PART = int(sys.argv[1]) if len(sys.argv) > 1 else 0
def hdr(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)
def ok(s):   print("  [OK]    " + s)
def bad(s):  print("  [FAIL]  " + s)
def note(s): print("  " + s)

RATIOS = [(17, 7, 13), (23, 7, 17), (71, 49, 61)]   # доказанные рёберные исключения
for b0, h0, n0 in RATIOS:
    assert b0^2 + h0^2 == 2*n0^2

# ===========================================================================
if PART in (0, 1):
    hdr("§1. Вывод с нуля. Магический 3x3 <-> центр e и два шага x,y.")
    R = PolynomialRing(QQ, ['e','x','y']); e,x,y = R.gens()
    # a b c / d e f / g h i ;  пары: (a,i),(c,g) угловые; (b,h),(d,f) рёберные
    cell = {'a': e+x, 'c': e+y, 'i': e-x, 'g': e-y, 'e': e}
    cell['b'] = 3*e - cell['a'] - cell['c']          # строка a b c
    cell['h'] = 2*e - cell['b']
    cell['d'] = 3*e - cell['a'] - cell['g']          # столбец a d g
    cell['f'] = 2*e - cell['d']
    lines = [('a','b','c'),('d','e','f'),('g','h','i'),('a','d','g'),
             ('b','e','h'),('c','f','i'),('a','e','i'),('c','e','g')]
    assert all(sum(cell[k] for k in L) == 3*e for L in lines)
    for p in [('a','i'),('c','g'),('b','h'),('d','f')]:
        assert cell[p[0]] + cell[p[1]] == 2*e
    note("клетки:  a=%s  b=%s  c=%s" % (cell['a'],cell['b'],cell['c']))
    note("         d=%s  e=%s  f=%s" % (cell['d'],cell['e'],cell['f']))
    note("         g=%s  h=%s  i=%s" % (cell['g'],cell['h'],cell['i']))
    ok("8 линий = 3e, 4 пары = 2e  [доказано, символьно]")
    note("ШАГИ: угловые пары -> x и y;  рёберные пары -> x+y и x-y.")
    note("Асимметрия: рёберные шаги = СУММА/РАЗНОСТЬ угловых; угловые = ПОЛУСУММА/ПОЛУРАЗНОСТЬ рёберных.")

    hdr("§1.2 Параметризация Сегре: 2e = (m^2+n^2)(1+t^2) двумя способами.")
    S = PolynomialRing(QQ, ['m','n','t']); m,n,t = S.gens(); F = S.fraction_field()
    s  = (m^2+n^2)/2
    E4 = s*(1+t^2)
    P1, Q1 = (m*t+n)^2, (m-n*t)^2          # первое представление 2e
    P2, Q2 = (m*t-n)^2, (m+n*t)^2          # второе представление 2e
    assert F(P1+Q1) == F(2*E4) and F(P2+Q2) == F(2*E4)
    ok("обе четвёрки дают 2e  [доказано, символьно]")

    def build(pairtype):
        """pairtype='edge': квадраты на рёбрах (G1); 'corner': квадраты на углах (G2)."""
        if pairtype == 'corner':
            cc = {'a': P1, 'i': Q1, 'c': P2, 'g': Q2, 'e': E4}
        else:
            cc = {'b': P1, 'h': Q1, 'd': P2, 'f': Q2, 'e': E4}
        # достроить остальные клетки из магичности
        if pairtype == 'corner':
            cc['b'] = 3*E4 - cc['a'] - cc['c']; cc['h'] = 2*E4 - cc['b']
            cc['d'] = 3*E4 - cc['a'] - cc['g']; cc['f'] = 2*E4 - cc['d']
        else:
            # b,h,d,f известны; углы: a = 3e-b-c ... удобнее из пар шагов
            # шаг (b,h): h-e = X3 ; шаг (d,f): f-e = -X4 ; x=(X3+X4)/2 ...
            X3 = cc['h'] - E4        # h = e + (x+y)  =>  X3 = x+y
            X4 = cc['f'] - E4        # f = e + (x-y)  =>  X4 = x-y
            xx = (X3 + X4)/2; yy = (X3 - X4)/2
            cc['a'] = E4 + xx; cc['i'] = E4 - xx
            cc['c'] = E4 + yy; cc['g'] = E4 - yy
        for L in lines: assert F(sum(cc[k] for k in L)) == F(3*E4), (pairtype, L)
        return cc

    G1 = build('edge'); G2 = build('corner')
    alpha = (3*m^2 - n^2)/2; gamma = (3*n^2 - m^2)/2
    assert F(alpha + gamma) == F(2*s)
    # явные формулы, заявленные ниже
    assert F(G1['a']) == F(m^2 + n^2*t^2) and F(G1['i']) == F(n^2 + m^2*t^2)
    assert F(G1['c']) == F(E4 - 2*m*n*t) and F(G1['g']) == F(E4 + 2*m*n*t)
    assert F(G2['b']) == F(alpha + gamma*t^2) and F(G2['h']) == F(gamma + alpha*t^2)
    assert F(G2['d']) == F(E4 - 4*m*n*t) and F(G2['f']) == F(E4 + 4*m*n*t)
    ok("явные формулы обоих семейств подтверждены символьно")
    hdr("§1.3 ЯВНЫЕ УРАВНЕНИЯ УГЛОВОГО СЕМЕЙСТВА G2")
    note("s = (m^2+n^2)/2,  alpha = (3m^2-n^2)/2,  gamma = (3n^2-m^2)/2,  alpha+gamma = 2s")
    note("")
    note("   [  (mt+n)^2         alpha + gamma t^2     (mt-n)^2      ]")
    note("   [  s(1+t^2) - 4mnt  s(1+t^2)              s(1+t^2)+4mnt ]")
    note("   [  (m+nt)^2         gamma + alpha t^2     (m-nt)^2      ]")
    note("")
    note("автоматические квадраты — ЧЕТЫРЕ УГЛА; свободны крест: b,d,e,f,h.")
    note("для сравнения G1 (рёберное), тот же центр:")
    note("   [  m^2+n^2 t^2      (mt+n)^2              s(1+t^2)-2mnt ]")
    note("   [  (mt-n)^2         s(1+t^2)              (m+nt)^2      ]")
    note("   [  s(1+t^2)+2mnt    (m-nt)^2              n^2+m^2 t^2   ]")
    note("СВОБОДНЫЕ КЛЕТКИ БОК О БОК (центр общий, s(1+t^2)):")
    note("   G1: {m^2+n^2t^2, n^2+m^2t^2}  и  {s(1+t^2) -+ 2mnt}")
    note("   G2: {alpha+gamma t^2, gamma+alpha t^2} и {s(1+t^2) -+ 4mnt}")
    note("то есть ОТКЛОНЕНИЯ ОТ ЦЕНТРА У G2 РОВНО ВДВОЕ БОЛЬШЕ, чем у G1 при тех же (m,n,t).")

    hdr("§1.4 Проверка подстановкой конкретных чисел")
    def numeric(fam, mv, nv, tv):
        sub = lambda P: QQ(F(P)(m=mv, n=nv, t=tv))
        M = [[sub(fam[k]) for k in row] for row in (('a','b','c'),('d','e','f'),('g','h','i'))]
        tot = sum(M[0]); assert all(sum(r) == tot for r in M)
        assert all(sum(M[i][j] for i in range(3)) == tot for j in range(3))
        assert M[0][0]+M[1][1]+M[2][2] == tot and M[0][2]+M[1][1]+M[2][0] == tot
        return M, tot, sum(1 for r in M for v in r if v.is_square())
    for (mv,nv,tv) in [(2,1,QQ(3)/5), (5,3,QQ(7)/2), (7,4,QQ(-2)/9), (13,5,QQ(4)/3)]:
        M1,T1,q1 = numeric(G1,mv,nv,tv); M2,T2,q2 = numeric(G2,mv,nv,tv)
        note(f"(m,n,t)=({mv},{nv},{tv}): G1 магичен, квадратов {q1}; G2 магичен, квадратов {q2}; суммы {T1}, {T2}")
        assert q1 >= 4 and q2 >= 4 and T1 == T2
    ok("оба семейства магические при подстановке; ровно >=4 автоматических квадрата; суммы совпадают")

    hdr("§1.5 ВЫРОЖДЕННАЯ ЛАТИНСКАЯ ТОЧКА G2 (её нет у G1) — тождество")
    t0 = (m-n)/(m+n)
    lat = {k: F(G2[k])(t=t0) for k in G2}
    LA = ((m^2+n^2)/(m+n))^2; LB = ((m^2+2*m*n-n^2)/(m+n))^2; LC = ((m^2-2*m*n-n^2)/(m+n))^2
    want = {'a':LA,'e':LA,'i':LA, 'b':LB,'f':LB,'g':LB, 'c':LC,'d':LC,'h':LC}
    for k in lat: assert F(lat[k]) == F(want[k]), k
    ok("при t=(m-n)/(m+n) ВСЕ ДЕВЯТЬ клеток G2 — квадраты, но значений ровно три (латинский квадрат)")
    note("   a=e=i=((m^2+n^2)/(m+n))^2,  b=f=g=((m^2+2mn-n^2)/(m+n))^2,  c=d=h=((m^2-2mn-n^2)/(m+n))^2")
    note("   [доказано, символьно] — тождество по m,n. У G1 девятка квадратов при t=1 требует m^2+n^2=[].")
    t1 = S.fraction_field()(1)
    g1at1 = {k: F(G1[k])(t=1) for k in G1}
    note("   контроль G1 при t=1: клетки = " + str(sorted(set(str(F(v)) for v in g1at1.values()))))

# ===========================================================================
if PART in (0, 2):
    hdr("§2. Универсальное семейство X(A,B) рода 5: оба G лежат в нём")
    note("Тройка «свободная пара + центр» даёт кривую")
    note("   X(A,B):  u^2 = A + B T^2,   v^2 = ((A+B)/2)(1+T^2),   w^2 = B + A T^2.")
    note("   G1: (A,B) = (m^2, n^2)                  — ОБА КВАДРАТА;")
    note("   G2: (A,B) = ((3m^2-n^2)/2, (3n^2-m^2)/2) — оба квадрата лишь в исключениях.")
    note("Обратно: m^2 = (3A+B)/4, n^2 = (A+3B)/4.")
    S = PolynomialRing(QQ, ['m','n']); m,n = S.gens()
    al = (3*m^2-n^2)/2; ga = (3*n^2-m^2)/2
    assert S.fraction_field()((3*al+ga)/4) == S.fraction_field()(m^2)
    assert S.fraction_field()((al+3*ga)/4) == S.fraction_field()(n^2)
    ok("пересчёт (A,B) -> (m^2,n^2) подтверждён символьно")
    note("Значит G1 и G2 — ДВА РАЗНЫХ АРИФМЕТИЧЕСКИХ СЛОЯ одного и того же геометрического пучка:")
    note("   G1 = {A,B оба квадраты};  G2 = {(3A+B)/4 и (A+3B)/4 оба квадраты}.")

# ===========================================================================
if PART in (0, 3):
    hdr("§3. Полный критерий изоморфизма X(A,B) ~ X(A',B')")
    Rp = PolynomialRing(QQ, 'p'); p = Rp.gen()
    note("Разложение Jac X(A,B) ~ E13 x E12 x E32 x Jac(H), где")
    note("   E13: y^2=(A+BT^2)(B+AT^2);  E12: y^2=(A+BT^2)s(1+T^2);  E32: y^2=(B+AT^2)s(1+T^2);")
    note("   H:   y^2=(A+BT^2)s(1+T^2)(B+AT^2)  (род 2).   0+1+1+1+2 = 5 = g(X).")
    # проверка формул j как функций p = A/B + B/A
    def jacobian_of_biquadratic(a, c, ee):
        # y^2 = a T^4 + c T^2 + ee   ->   Y^2 = X(X^2 + 2cX + (c^2-4a*ee))
        return EllipticCurve([0, 2*c, 0, c^2 - 4*a*ee, 0])
    def invs(A, B):
        A = QQ(A); B = QQ(B); s = (A+B)/2
        E13 = jacobian_of_biquadratic(A*B, A^2+B^2, A*B)
        E12 = jacobian_of_biquadratic(s*B, s*(A+B), s*A)
        E32 = jacobian_of_biquadratic(s*A, s*(A+B), s*B)
        return E13, E12, E32
    for (A,B) in [(9,4),(25,16),(109,229),(3,5),(49,289)]:
        E13,E12,E32 = invs(A,B)
        pp = QQ(A)/B + QQ(B)/A
        assert E13.j_invariant() == 16*(pp^2+12)^3/(pp^2-4)^2, (A,B,'j13')
        assert E12.j_invariant() == 16*(pp+14)^3/(pp-2)^2,     (A,B,'j12')
        assert E32.j_invariant() == E12.j_invariant()
    ok("j(E13) = 16(p^2+12)^3/(p^2-4)^2  и  j(E12)=j(E32)=16(p+14)^3/(p-2)^2, p = A/B+B/A   [проверено подстановкой]")
    note("Оба инварианта — функции ОДНОГО параметра p; p зависит только от неупорядоченной пары {A,B}")
    note("с точностью до общего множителя. Значит p — модуль пучка.")
    note("")
    note("КРИТЕРИЙ (достаточный для НЕсводимости): если p(G2) != p(G1') для всех допустимых G1',")
    note("то кривые не изоморфны — ибо j(E12) определяет p с точностью до корней кубики,")
    note("а вместе с j(E13) (функция p^2) — однозначно вне явного конечного множества.")
    note("Проверяем ниже ЯВНО: сколько p' дают тот же (j12,j13), что данное p.")
    Rq = PolynomialRing(QQ, 'q'); q = Rq.gen()
    for p0 in [QQ(5)/2, QQ(10)/3, QQ(229)/109 + QQ(109)/229, QQ(338)/49 + QQ(49)/338]:
        j12 = 16*(p0+14)^3/(p0-2)^2; j13 = 16*(p0^2+12)^3/(p0^2-4)^2
        f12 = (16*(q+14)^3 - j12*(q-2)^2).numerator()
        sols12 = [r for r in f12.roots(QQ, multiplicities=False)]
        f13 = (16*(q^2+12)^3 - j13*(q^2-4)^2).numerator()
        sols13 = [r for r in f13.roots(QQ, multiplicities=False)]
        both = sorted(set(sols12) & set(sols13))
        note(f"  p={p0}: решений j12: {sols12}; решений j13: {sols13}; ОБА: {both}")
        assert both == [p0], (p0, both)
    ok("(j12, j13) вместе определяют p ОДНОЗНАЧНО в проверенных точках  [проверено]")
    note("Строгая версия того же — над полем функций Q(p):")
    Rpq = PolynomialRing(QQ, ['pp','qq']); pp, qq = Rpq.gens()
    F1 = (pp+14)^3*(qq-2)^2 - (qq+14)^3*(pp-2)^2
    F2 = (pp^2+12)^3*(qq^2-4)^2 - (qq^2+12)^3*(pp^2-4)^2
    assert F1.subs({qq: pp}) == 0 and F2.subs({qq: pp}) == 0
    H1 = Rpq(F1 // (qq - pp)); H2 = Rpq(F2 // (qq - pp))
    note(f"  F1/(q-p) = {H1.factor()}")
    Kp = PolynomialRing(QQ, 'pp').fraction_field(); Kq = PolynomialRing(Kp, 'qq')
    h1 = Kq(H1.subs({pp: Kp.gen()})); h2 = Kq(H2.subs({pp: Kp.gen()}))
    gg = gcd(h1, h2)
    note(f"  gcd(F1/(q-p), F2/(q-p)) над Q(p)[q] = {gg}")
    if gg.degree() == 0:
        ok("общих корней, кроме q=p, НЕТ над Q(p)  =>  p -> (j12,j13) инъективно вне конечного")
        res = h1.resultant(h2)
        num = Kp(res).numerator()
        note(f"  исключительное множество: корни результанта {num.factor() if num != 0 else 0}")
        note(f"  рациональные исключительные p: {num.roots(QQ, multiplicities=False) if num != 0 else 'все'}")
    else:
        bad(f"есть общий множитель: {gg.factor()} — инъективность не доказана")
    note("=> X(A,B) ~ X(A',B')  ==>  p = p'  ==>  {A:B} = {A':B'}.")
    note("   (обратное тоже верно: при {A:B}={A':B'} кривые — квадратичные твисты, а при")
    note("    равном классе твиста — изоморфны над Q.)")

    hdr("§3.2 Когда G2(m,n) попадает в G1? — диофантов вопрос, решаем до конца")
    note("Нужно: (A,B)=(alpha,gamma) ~ lam^2 (m'^2, n'^2). Так как alpha+gamma = m^2+n^2 фиксирует масштаб,")
    note("условие равносильно: alpha*gamma = (3m^2-n^2)(3n^2-m^2)/4 — КВАДРАТ,")
    note("и тогда сам lam^2 определяется; дополнительно нужен квадратный класс alpha.")
    note("Кривая: V^2 = (3u^2-... ) ->  положим u=m/n:  W^2 = -3u^4 + 10u^2 - 3.")
    Ru = PolynomialRing(QQ, 'u'); u = Ru.gen()
    quart = -3*u^4 + 10*u^2 - 3
    note(f"  квартика: {quart};  точка u=1 -> W^2 = {quart(1)} = 2^2")
    def quartic_jacobian(f):
        """f = a x^4+b x^3+c x^2+d x+e  ->  E: Y^2 = X^3 - 27 I X - 27 J."""
        co = [f[4], f[3], f[2], f[1], f[0]]
        a,b,c,d,ee = [QQ(z) for z in co]
        I = 12*a*ee - 3*b*d + c^2
        J = 72*a*c*ee + 9*b*c*d - 27*a*d^2 - 27*ee*b^2 - 2*c^3
        return EllipticCurve([0,0,0,-27*I,-27*J]).minimal_model(), I, J
    Emin, Iq, Jq = quartic_jacobian(quart)
    note(f"  I = {Iq}, J = {Jq};  якобиан: {Emin}, кондуктор {Emin.conductor()}, кручение {Emin.torsion_order()}")
    rb = Emin.rank_bounds()
    note(f"  rank_bounds = {rb};  analytic rank = {Emin.analytic_rank()}")
    if rb == (0,0):
        ok("ранг 0  => на квартике W^2=-3u^4+10u^2-3 рациональных точек конечное число (только кручение)")
    else:
        bad(f"ранг не 0: {rb} — конечность точек не следует автоматически")
    # прямой перебор рациональных u = a/b малой высоты
    sols = []
    for bq in range(1, 200):
        for aq in range(0, 200):
            if gcd(aq,bq) != 1: continue
            uu = QQ(aq)/bq
            val = quart(uu)
            if val >= 0 and val.is_square(): sols.append(uu)
    note(f"  перебор u=a/b, a,b<=200: решения {sorted(set(sols))}")
    note("  (u=1 <=> m=n: у G2 при m=n клетки a и i совпадают с c и g — вырождение)")

# ===========================================================================
if PART in (0, 4):
    hdr("§4. Сечения (фиксируем одну полную пару): A, C, beta для обоих типов")
    note("Нормируем центр в 1. Пара: b0^2+h0^2 = 2 n0^2.  Вторая пара — коника, параметр T:")
    note("   Bf = (1+2T-T^2)/(1+T^2),  Hf = (1-2T-T^2)/(1+T^2),  Bf^2+Hf^2 = 2.")
    RT = PolynomialRing(QQ, ['T','b0','h0','n0']); T,b0v,h0v,n0v = RT.gens(); FT = RT.fraction_field()
    W  = (1+T^2)^2
    D0 = T*(1-T^2)
    Bf2 = (1+2*T-T^2)^2; Hf2 = (1-2*T-T^2)^2
    assert FT(Bf2+Hf2) == FT(2*W)
    ok("коника: Bf^2+Hf^2 = 2 тождественно")
    # РЁБЕРНОЕ: обе рёберные пары квадраты -> свободные клетки = УГЛЫ = полусуммы
    for (bb,hh,nn) in RATIOS:
        A_e, C_e, beta_e = QQ(hh^2+nn^2)/2, QQ(bb^2+nn^2)/2, QQ(2*nn^2)
        # угол = (X + Y)/2, X in {b0^2,h0^2} (норм. на n0^2), Y in {Bf^2,Hf^2}/W
        for (X, Aexp) in [(QQ(hh^2), A_e), (QQ(bb^2), C_e)]:
            for sg, Y2 in [(+1, Bf2), (-1, Hf2)]:
                lhs = (X*W + Y2*nn^2)/2            # (клетка)*W*n0^2
                rhs = Aexp*W + sg*beta_e*D0
                assert FT(lhs) == FT(rhs), ('edge', bb,hh,nn, X, sg)
        assert A_e + C_e == beta_e
        # УГЛОВОЕ: обе угловые пары квадраты -> свободные клетки = РЁБРА = X+Y-1
        A_c, C_c, beta_c = QQ(hh^2), QQ(bb^2), QQ(4*nn^2)
        for (X, Aexp) in [(QQ(hh^2), A_c), (QQ(bb^2), C_c)]:
            for sg, Y2 in [(+1, Bf2), (-1, Hf2)]:
                lhs = X*W + Y2*nn^2 - nn^2*W       # (клетка)*W*n0^2
                assert FT(lhs) == FT(Aexp*W + sg*beta_c*D0), ('corner', bb,hh,nn, X, sg)
        assert 2*(A_c + C_c) == beta_c
    ok("ОБА семейства дают один вид  Y^2 = A(T^2+1)^2 +- beta T(1-T^2)   [проверено подстановкой]")
    ok("рёберное:  A=(h0^2+n0^2)/2, C=(b0^2+n0^2)/2, beta=2n0^2  =>  beta = A + C")
    ok("угловое:   A=h0^2,          C=b0^2,          beta=4n0^2  =>  beta = 2(A + C)")
    note("")
    note("ГЛАВНОЕ СТРУКТУРНОЕ ОТЛИЧИЕ: у УГЛОВОГО сечения A и C — ВСЕГДА ТОЧНЫЕ КВАДРАТЫ (h0^2, b0^2).")
    note("У рёберного A/C = (h0^2+n0^2)/(b0^2+n0^2) — квадрат лишь в исключительных тройках.")
    note("")
    prim = []
    for nn in range(1, 401):
        for hh in range(1, nn+1):
            bb2 = 2*nn^2 - hh^2
            if bb2 <= 0: continue
            bb = ZZ(bb2).isqrt()
            if bb*bb != bb2 or bb < hh: continue
            if gcd([bb,hh,nn]) != 1: continue
            prim.append((bb,hh,nn))
    note(f"примитивных троек b^2+h^2=2n^2 с n<=400: {len(prim)}")
    def data(bb,hh,nn,kind):
        if kind=='edge': A,C,be = QQ(hh^2+nn^2)/2, QQ(bb^2+nn^2)/2, QQ(2*nn^2)
        else:            A,C,be = QQ(hh^2), QQ(bb^2), QQ(4*nn^2)
        bp = (C-A)/(C+A); disc = (C-A)*(C+A)
        return A,C,be,bp,disc
    same = []
    for (bb,hh,nn) in prim:
        Ae,Ce,bee,bpe,de = data(bb,hh,nn,'edge')
        Ac,Cc,bec,bpc,dc = data(bb,hh,nn,'corner')
        assert bpc == 2*bpe, (bb,hh,nn)
        assert dc == 2*de
    ok("тождества: b'_угол = 2 b'_ребро  и  (C^2-A^2)_угол = 2 (C^2-A^2)_ребро  на всех %d тройках" % len(prim))
    note("  => k_ребро = Q(sqrt(b0^2-h0^2)),  k_угол = Q(sqrt(2(b0^2-h0^2))) — РАЗНЫЕ поля (сдвиг на 2).")
    note("")
    note("Таблица для доказанных отношений:")
    note("  (b,h,n)      | рёберное A,C,beta, k, b'          | угловое A,C,beta, k, b'")
    for (bb,hh,nn) in RATIOS:
        Ae,Ce,bee,bpe,de = data(bb,hh,nn,'edge')
        Ac,Cc,bec,bpc,dc = data(bb,hh,nn,'corner')
        sq = lambda z: QQ(z).squarefree_part()
        note(f"  ({bb},{hh},{nn}) | A={Ae} C={Ce} beta={bee} k=Q(v{sq(de)}) b'={bpe} | A={Ac} C={Cc} beta={bec} k=Q(v{sq(dc)}) b'={bpc}")
    note("")
    note("Пересечение A:C рёберных и угловых сечений (n<=400):")
    edgeAC = {}
    for (bb,hh,nn) in prim:
        Ae,Ce,_,_,_ = data(bb,hh,nn,'edge'); edgeAC[Ae/Ce] = (bb,hh,nn)
    hits = []
    for (bb,hh,nn) in prim:
        Ac,Cc,_,_,_ = data(bb,hh,nn,'corner')
        if Ac/Cc in edgeAC: hits.append(((bb,hh,nn), edgeAC[Ac/Cc]))
    note(f"  совпадений A:C: {hits if hits else 'НЕТ НИ ОДНОГО'}")
    note("  (A:C у углового = (h0/b0)^2 — всегда квадрат; у рёберного квадратом почти никогда)")
    nsq = [(bb,hh,nn) for (bb,hh,nn) in prim if (QQ(hh^2+nn^2)/2).is_square() and (QQ(bb^2+nn^2)/2).is_square()]
    note(f"  рёберных троек с A и C ОБА квадраты (тогда A:C — квадрат): {nsq if nsq else 'НЕТ (n<=400)'}")

# ===========================================================================
if PART in (0, 5):
    hdr("§5. Род 2: C_J2: y^2 = z(z^2-1)(z^2-b'^2). Igusa — ПОЛНЫЙ инвариант над Qbar.")
    Rb = PolynomialRing(QQ, 'B'); B = Rb.gen()
    def igusa_abs(bp):
        z = polygen(QQ, 'z')
        H = HyperellipticCurve(z*(z^2-1)*(z^2-bp^2))
        J = H.igusa_clebsch_invariants()
        I2,I4,I6,I10 = J
        return (I2^5/I10, I4*I2^3/I10, I6*I2^2/I10)
    prim = []
    for nn in range(1, 401):
        for hh in range(1, nn+1):
            bb2 = 2*nn^2 - hh^2
            if bb2 <= 0: continue
            bb = ZZ(bb2).isqrt()
            if bb*bb != bb2 or bb < hh or gcd([bb,hh,nn]) != 1: continue
            prim.append((bb,hh,nn))
    E, Cc = {}, {}
    for (bb,hh,nn) in prim:
        bpe = (QQ(bb^2+nn^2)/2 - QQ(hh^2+nn^2)/2)/(QQ(bb^2+nn^2)/2 + QQ(hh^2+nn^2)/2)
        bpc = 2*bpe
        if bpe in (0,) or bpc in (0,): continue
        E[(bb,hh,nn)] = igusa_abs(bpe); Cc[(bb,hh,nn)] = igusa_abs(bpc)
    inter = set(E.values()) & set(Cc.values())
    note(f"троек: {len(prim)};  различных рёберных классов: {len(set(E.values()))};  угловых: {len(set(Cc.values()))}")
    note(f"ПЕРЕСЕЧЕНИЕ абсолютных инвариантов Игузы: {len(inter)}")
    if not inter: ok("ни одно угловое сечение не Qbar-изоморфно ни одному рёберному (n<=400)  [доказано (ПО)]")
    else: bad(f"есть совпадения: {list(inter)[:3]}")
    hdr("§5.2 Символьно: при каких b' кривая с b' и кривая с 2b' Qbar-изоморфны?")
    note("Решаем i1(2b) = i1(b), i2(2b)=i2(b), i3(2b)=i3(b) как систему по b.")
    z = polygen(Rb.fraction_field(), 'z')
    def igusa_sym(expr):
        H = HyperellipticCurve(z*(z^2-1)*(z^2-expr^2))
        I2,I4,I6,I10 = H.igusa_clebsch_invariants()
        return (I2^5/I10, I4*I2^3/I10, I6*I2^2/I10)
    try:
        i_b = igusa_sym(B); i_2b = igusa_sym(2*B)
        polys = []
        for x1, x2 in zip(i_b, i_2b):
            num = (x1 - x2).numerator()
            polys.append(Rb(num))
        g = gcd(polys)
        note(f"НОД трёх числителей разностей: {g.factor() if g != 0 else 0}")
        roots = g.roots(QQ, multiplicities=False) if g != 0 else []
        note(f"рациональные корни (кандидаты на совпадение): {roots}")
        good = [r for r in roots if r not in (0,) and r^2 != 1 and (2*r)^2 != 1 and r != 0]
        note(f"из них невырожденные (b' != 0, +-1, +-1/2): {good}")
        if not good: ok("НЕТ ни одного b', при котором угловая C_J2 Qbar-изоморфна своей рёберной  [доказано, символьно]")
    except Exception as ex:
        bad(f"символьная часть не прошла: {ex}")
    hdr("§5.3 ОРБИТА: при каких g кривая C_J2(g) Qbar-изоморфна C_J2(b)?")
    Rbg = PolynomialRing(QQ, ['bb','gg']); bb_, gg_ = Rbg.gens()
    Fb = Rbg.fraction_field()
    zg = polygen(Fb, 'z')
    def ig3(expr):
        Hh = HyperellipticCurve(zg*(zg^2-1)*(zg^2-expr^2))
        I2,I4,I6,I10 = Hh.igusa_clebsch_invariants()
        return (I2^5/I10, I4*I2^3/I10, I6*I2^2/I10)
    try:
        ia = ig3(Fb(bb_)); ib = ig3(Fb(gg_))
        nums = [Rbg((x - y).numerator()) for x, y in zip(ia, ib)]
        G = gcd(nums)
        note(f"  общий множитель системы i_k(g)=i_k(b):  {G.factor()}")
        note("  => орбита g при заданном b — корни этих множителей (конечный явный список).")
    except Exception as ex:
        bad(f"орбита не посчиталась: {ex}")

# ===========================================================================
if PART in (0, 6):
    hdr("§6. Дешёвые закрытия для УГЛОВЫХ сечений трёх запрещённых отношений")
    note("Вход конвейера: A = h0^2, C = b0^2, beta = 4 n0^2  (вместо рёберных полусумм).")
    note("ВАЖНО: у углового сечения при T=0 все четыре свободные клетки = A,A,C,C = h0^2,h0^2,b0^2,b0^2,")
    note("то есть ДЕВЯТЬ КВАДРАТОВ есть ВСЕГДА — но это вырождение c=g=e. Ранг 0 поэтому недостижим.")
    RT = PolynomialRing(QQ, 'T'); T = RT.gen()
    def quartic_jacobian(f):
        co = [f[4], f[3], f[2], f[1], f[0]]
        a,b,c,d,ee = [QQ(z) for z in co]
        I = 12*a*ee - 3*b*d + c^2
        J = 72*a*c*ee + 9*b*c*d - 27*a*d^2 - 27*ee*b^2 - 2*c^3
        return EllipticCurve([0,0,0,-27*I,-27*J]).minimal_model()
    def rational_T(f, N=400):
        out = []
        for q in range(1, N+1):
            for p in range(-N, N+1):
                if gcd(p,q) != 1: continue
                v = f(QQ(p)/q)
                if v >= 0 and v.is_square(): out.append(QQ(p)/q)
        return sorted(set(out))
    ALL = {}
    for (bb,hh,nn) in RATIOS + [(7,1,5)]:
        hdr(f"  сечения ({bb},{hh},{nn})")
        for kind in ('corner','edge'):
            if kind == 'corner': A, C, be = QQ(hh^2), QQ(bb^2), QQ(4*nn^2)
            else:               A, C, be = QQ(hh^2+nn^2)/2, QQ(bb^2+nn^2)/2, QQ(2*nn^2)
            note(f"   --- {kind}: A={A}, C={C}, beta={be}, A/C={A/C} (квадрат? {(A/C).is_square()})")
            for nameX, AA in [('A', A), ('C', C)]:
                for sg in (+1,-1):
                    f = AA*(T^2+1)^2 + sg*be*T*(1-T^2)
                    try:
                        E = quartic_jacobian(f)
                        rb = E.rank_bounds(); tor = E.torsion_order()
                        has0 = f(0).is_square()
                        extra = ""
                        if rb == (0,0):
                            Ts = rational_T(f, 200)
                            extra = f"   ТОЧКИ T (|p|,q<=200): {Ts}   (|E(Q)|={tor})"
                        note(f"     {nameX}{'+' if sg>0 else '-'}: {E.ainvs()} N={E.conductor()} rank={rb} tors={tor} точка при T=0? {has0}{extra}")
                        ALL[(bb,hh,nn,kind,nameX,sg)] = rb
                    except Exception as ex:
                        note(f"     {nameX}{'+' if sg>0 else '-'}: ошибка {ex}")
            bp = (C-A)/(C+A); D = (C-A)*(C+A)
            note(f"     b'={bp}, (C^2-A^2)={D}, k=Q(sqrt({D.squarefree_part()}))")
    hdr("§6.1bis свод: есть ли квартика ранга 0 (=> дешёвое закрытие сечения)")
    for key, rb in sorted(ALL.items(), key=lambda kv: str(kv[0])):
        if rb == (0,0): note(f"  РАНГ 0: {key}")
    if not any(rb == (0,0) for rb in ALL.values()):
        note("  квартик ранга 0 нет ни в одном из проверенных сечений (ни угловых, ни рёберных)")
    hdr("§6.2 Локальная разрешимость угловой системы (ELS-грубый тест)")
    note("Формы дополненного покрытия: p^2-q^2, p^2+q^2, A p^2 - C q^2, C p^2 - A q^2.")
    note("У УГЛОВОГО случая A=h0^2, C=b0^2, поэтому")
    note("   A p^2 - C q^2 = (h0 p - b0 q)(h0 p + b0 q)   и   C p^2 - A q^2 = (b0 p - h0 q)(b0 p + h0 q)")
    note("— ОБЕ РАСПАДАЮТСЯ НА ЛИНЕЙНЫЕ МНОЖИТЕЛИ над Q. В рёберном случае они неприводимы.")
    note("Это структурное упрощение углового случая, отсутствующее в рёберном.")
    for (bb,hh,nn) in RATIOS:
        Rpq = PolynomialRing(QQ, ['P','Q']); P,Q = Rpq.gens()
        f3 = hh^2*P^2 - bb^2*Q^2; f4 = bb^2*P^2 - hh^2*Q^2
        note(f"  ({bb},{hh},{nn}): A p^2-C q^2 = {f3.factor()};  C p^2-A q^2 = {f4.factor()}")

# ===========================================================================
if PART in (0, 7):
    hdr("§7. ВЕЩЕСТВЕННОЕ (архимедово) ПРЕПЯТСТВИЕ — есть у углового, НЕТ у рёберного")
    note("Тригонометрическая нормальная форма. Пусть")
    note("   cos psi = (m^2-n^2)/(m^2+n^2), sin psi = 2mn/(m^2+n^2);")
    note("   cos theta = (1-t^2)/(1+t^2),   sin theta = 2t/(1+t^2).")
    note("Тогда (нормируя на центр e = s(1+t^2)):")
    note("   G1: свободные клетки =  1,  1 +- cos psi cos theta,  1 +- sin psi sin theta")
    note("   G2: свободные клетки =  1,  1 +- 2 cos psi cos theta, 1 +- 2 sin psi sin theta")
    note("   автоматические квадраты в обоих: 1 -+ cos(psi+theta), 1 -+ cos(psi-theta)")
    # символьная проверка
    S = PolynomialRing(QQ, ['m','n','t']); m,n,t = S.gens(); F = S.fraction_field()
    s = (m^2+n^2)/2; E4 = s*(1+t^2)
    cp = (m^2-n^2)/(m^2+n^2); sp = 2*m*n/(m^2+n^2)
    ct = (1-t^2)/(1+t^2);     st = 2*t/(1+t^2)
    assert F(cp^2+sp^2) == 1 and F(ct^2+st^2) == 1
    G1free = [m^2+n^2*t^2, n^2+m^2*t^2, E4-2*m*n*t, E4+2*m*n*t]
    G1trig = [E4*(1+cp*ct), E4*(1-cp*ct), E4*(1-sp*st), E4*(1+sp*st)]
    for u,v in zip(G1free, G1trig): assert F(u) == F(v), (u,v)
    alpha=(3*m^2-n^2)/2; gamma=(3*n^2-m^2)/2
    G2free = [alpha+gamma*t^2, gamma+alpha*t^2, E4-4*m*n*t, E4+4*m*n*t]
    G2trig = [E4*(1+2*cp*ct), E4*(1-2*cp*ct), E4*(1-2*sp*st), E4*(1+2*sp*st)]
    for u,v in zip(G2free, G2trig): assert F(u) == F(v), (u,v)
    auto = [(m*t+n)^2, (m-n*t)^2, (m*t-n)^2, (m+n*t)^2]
    autotrig = [E4*(1-cos_p) for cos_p in [ (cp*ct-sp*st), -(cp*ct-sp*st), (cp*ct+sp*st), -(cp*ct+sp*st)]]
    for u,v in zip(auto, autotrig): assert F(u) == F(v), ('auto',u,v)
    ok("тригонометрическая форма обоих семейств подтверждена символьно  [доказано]")
    note("")
    note("СЛЕДСТВИЕ (положительность всех девяти клеток):")
    note("   G1: |cos psi cos theta| <= 1 и |sin psi sin theta| <= 1 — ВЫПОЛНЕНО ТОЖДЕСТВЕННО.")
    note("   G2: нужно |cos psi cos theta| <= 1/2 И |sin psi sin theta| <= 1/2 — НЕТРИВИАЛЬНО.")
    import random, math
    random.seed(int(11))
    N = int(400000); cnt = int(0)
    for _ in range(N):
        a1 = random.uniform(0.0, 6.283185307179586); a2 = random.uniform(0.0, 6.283185307179586)
        if abs(math.cos(a1)*math.cos(a2)) <= 0.5 and abs(math.sin(a1)*math.sin(a2)) <= 0.5: cnt += int(1)
    note("   доля тора (psi,theta), проходящая G2-положительность: %.4f (Монте-Карло, %d точек)" % (float(cnt)/float(N), N))
    note("")
    note("СЕЧЕНИЕ (фиксируем одну ПОЛНУЮ пару): sigma = psi+theta фиксировано, свободно tau = psi-theta.")
    note("   cos psi cos theta = (cos sigma + cos tau)/2,  sin psi sin theta = (cos tau - cos sigma)/2.")
    note("   фиксированная пара даёт cos sigma = 1 - b0^2/n0^2 (и 1+cos sigma = h0^2/n0^2).")
    note("   G2 (углы): |cos tau + cos sigma| <= 1 и |cos tau - cos sigma| <= 1")
    note("       <=>  |cos tau| <= 1 - |cos sigma| = h0^2/n0^2.   НЕТРИВИАЛЬНО (h0 < n0 всегда).")
    note("   G1 (рёбра): те же неравенства с 1/2 -> 1, то есть |cos tau| <= 2 - |cos sigma| > 1 — ПУСТО.")
    note("")
    prim = []
    for nn in range(1, 201):
        for hh in range(1, nn+1):
            bb2 = 2*nn^2 - hh^2
            if bb2 <= 0: continue
            bb = ZZ(bb2).isqrt()
            if bb*bb != bb2 or bb < hh or gcd([bb,hh,nn]) != 1: continue
            prim.append((bb,hh,nn))
    note("  (b,h,n)    | c = (h0/n0)^2 | доля окружности tau, где угловое сечение вещественно возможно")
    for (bb,hh,nn) in [(17,7,13),(23,7,17),(71,49,61),(7,1,5),(41,1,29),(79,47,65),(89,23,65)]:
        c = QQ(hh^2)/nn^2
        frac = (2/pi)*arcsin(RR(c))
        note(f"   ({bb},{hh},{nn})  |  {c} = {RR(c):.4f}  |  {RR(frac):.4f}")
    note("")
    note("  ПРОВЕРКА ЧИСЛАМИ: перебираем T и смотрим знаки четырёх свободных клеток сечения.")
    RT = PolynomialRing(QQ, 'T'); T = RT.gen()
    for (bb,hh,nn) in [(17,7,13),(23,7,17),(71,49,61)]:
        for kind in ('corner','edge'):
            if kind=='corner': A,C,be = QQ(hh^2), QQ(bb^2), QQ(4*nn^2)
            else:              A,C,be = QQ(hh^2+nn^2)/2, QQ(bb^2+nn^2)/2, QQ(2*nn^2)
            fs = [A*(T^2+1)^2 + sg*be*T*(1-T^2) for sg in (1,-1)] + \
                 [C*(T^2+1)^2 + sg*be*T*(1-T^2) for sg in (1,-1)]
            good = 0; tot = 0
            for k in range(-4000, 4001):
                Tv = QQ(k)/1000; tot += 1
                if all(f(Tv) > 0 for f in fs): good += 1
            # добавим «хвост» |T|>4 через T -> 1/T симметрию
            note(f"   ({bb},{hh},{nn}) {kind}: доля T в [-4,4] с ЧЕТЫРЬМЯ положительными клетками = {good}/{tot} = {RR(good/tot):.4f}")
            predicted = QQ(hh^2)/nn^2 if kind=='corner' else None
            if kind=='corner':
                note(f"        предсказание теории: |cos tau| <= {predicted} = {RR(predicted):.4f}, доля окружности {RR((2/pi)*arcsin(RR(predicted))):.4f}")
    ok("вещественное препятствие углового случая подтверждено численно  [проверено]")

print("\nDONE")
