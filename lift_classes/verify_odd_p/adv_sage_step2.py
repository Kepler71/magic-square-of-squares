#!/usr/bin/env python3
"""Независимая (другой движок: Sage Rational.is_padic_square) адверсариальная проверка Шага 2
теоремы Codex об остаточных классах.  Запуск:
    env DOT_SAGE=/tmp/claude_oddp_sage python3 adv_sage_step2.py 3,5,7,13 [TIME_LIMIT_SEC]

Для каждой (p; r,s) строим z, НАМЕРЕННО попадающие в опасные места:
  * z = -1/lam + t p^k  (lam из 8 ненулевых коэффициентов, k=-3..9):  клетка 1+lam z имеет
    оценку v(lam)+k+v(t) — кратные нули и в целом, и в полюсном режиме (p | lam => -1/lam полюс);
  * пары (r,s) с p^a | r, p^a | s, p^a | (r-s), p^a | (r+s), a=1..5, плюс малые пары;
  * случайные z = u p^e, e=-5..3.
Если все 8 произведений строк/столбцов/диагоналей арифметической сетки — квадраты в Q_p,
проверяются утверждения Шага 1-2 (и локальные следствия Шагов 4-5 при p=3, p=3 mod 4):
  C1 противоположные клетки: произведение — квадрат в Q_p (одинаковый полный класс);
  C2 вектор чётностей = матрица (S,U,S+U);
  C3 z in Z_p => все чётности 0;
  C4 v(T) нечётна => p|(r-s);  v(L) нечётна => p|(r+s)   (T,L считаются НАПРЯМУЮ как произведения);
  C5 p|r или p|s, полюс => n чётно, S=U=0;
  C6 полюс, r,s единицы: n нечётно => U=1 и ровно одна из r+-s кратна p и S=[p|(r-s)]; n чётно => S=U=0;
  C7 q_lam = (n-v(lam)) mod 2 при v(lam)<n, иначе 0;
  C8 p=3 mod 4 => все чётности 0;  p=3 => v_3(z)>=1.
Отрицательные контроли (должны НАРУШАТЬСЯ, иначе проверка пуста):
  N1 перепутанная связь «v(T) нечётна => p|(r+s)»;
  N2 без двух краевых линий (только 6 условий) C4 уже не обязано выполняться.
Все циклы конечны; есть лимит времени.
"""
import sys, time, random, itertools
from sage.all import QQ, ZZ, gcd

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
LINES = ([[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] +
         [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] +
         [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
EDGE = [LINES[0], LINES[2]]          # строки i=-1, i=+1 (для контроля N2 выбрасываем их)
PAIRS = {'r': ((1, 0), (-1, 0)), 's': ((0, 1), (0, -1)),
         'p': ((1, 1), (-1, -1)), 'm': ((1, -1), (-1, 1))}


def mat(S, U):
    return {(-1, -1): S, (-1, 0): U, (-1, 1): (S + U) % 2,
            (0, -1): U, (0, 0): 0, (0, 1): U,
            (1, -1): (S + U) % 2, (1, 0): U, (1, 1): S}


MATS = [mat(S, U) for S in (0, 1) for U in (0, 1)]


def vp(x, p):
    return ZZ(x).valuation(p) if x != 0 else 10**9


def pairs_for(p, rng):
    out = set()
    for r in range(-9, 10):
        for s in range(1, 10):
            out.add((r, s))
    for a in range(1, 6):
        q = p ** a
        for m in (1, 2, 3, 4, 7, 11):
            for t in (1, 2, 3, 5, 8, 13):
                for sg in (1, -1):
                    Q = sg * q * m
                    out.update([(Q, t), (t, Q), (t + Q, t), (Q - t, t), (t, t + Q), (t, Q - t)])
    for _ in range(60):
        out.add((rng.randint(-10**4, 10**4), rng.randint(1, 10**4)))
    good = []
    for (r, s) in out:
        if r == 0 or s == 0 or r == s or r == -s or gcd(r, s) != 1:
            continue
        good.append((r, s))
    return sorted(good)


def zs_for(p, r, s, rng):
    lams = [r, -r, s, -s, r + s, -r - s, r - s, s - r]
    out = []
    for lam in lams:
        z0 = QQ(-1) / lam
        for k in range(-3, 10):
            for t in (1, -1, 2, 3, p + 1, 1 + 2 * p):
                out.append(z0 + QQ(t) * QQ(p) ** k)
    for e in range(-5, 4):
        for _ in range(4):
            u = rng.randint(1, p ** 7)
            if u % p == 0:
                u += 1
            out.append(QQ(u) * QQ(p) ** e * rng.choice((1, -1)))
    return out


def run(p, tlimit):
    rng = random.Random(20260926 + p)
    t0 = time.time()
    st = dict(tested=0, passed=0, pole=0, S1=0, SU1=0, multizero_int=0, multizero_pole=0,
              eq_n=0, pr_ps_pole=0, N1_viol=0, N2_pass=0, N2_viol=0)
    prs = pairs_for(p, rng)
    for idx, (r, s) in enumerate(prs):
        if time.time() - t0 > tlimit:
            print("  p=%d TIME LIMIT at pair %d/%d" % (p, idx, len(prs)), flush=True)
            break
        vr, vs, vpl, vmi = vp(r, p), vp(s, p), vp(r + s, p), vp(r - s, p)
        vl = {'r': vr, 's': vs, 'p': vpl, 'm': vmi}
        for z in zs_for(p, r, s, rng):
            if z == 0:
                continue
            f = {c: 1 + (c[0] * r + c[1] * s) * z for c in CELLS}
            if any(f[c] == 0 for c in CELLS):
                continue
            st['tested'] += 1
            sq = [QQ(f[a] * f[b] * f[c]).is_padic_square(p) for (a, b, c) in LINES]
            par = {c: f[c].valuation(p) % 2 for c in CELLS}
            T = f[(-1, 0)] * f[(1, 1)] * f[(0, -1)]
            L = f[(-1, 0)] * f[(1, -1)] * f[(0, 1)]
            vT, vL = T.valuation(p) % 2, L.valuation(p) % 2
            # N2: только 6 линий (без строк i=+-1)
            if all(sq[k] for k in range(8) if k not in (0, 2)):
                st['N2_pass'] += 1
                if (vT and (r - s) % p) or (vL and (r + s) % p):
                    st['N2_viol'] += 1
            if not all(sq):
                continue
            st['passed'] += 1
            ez = z.valuation(p)
            n = -ez
            # C1
            for c in CELLS:
                assert QQ(f[c] * f[(-c[0], -c[1])]).is_padic_square(p), ("C1", p, r, s, z, c)
            # C2
            assert any(all(par[c] == M[c] for c in CELLS) for M in MATS), ("C2", p, r, s, z, par)
            S, U, SU = par[(1, 1)], par[(1, 0)], par[(1, -1)]
            assert vT == S and vL == SU, ("TL", p, r, s, z)
            vals = [f[c].valuation(p) for c in CELLS]
            # C3
            if ez >= 0:
                if max(vals) >= 2:
                    st['multizero_int'] += 1
                assert all(x == 0 for x in par.values()), ("C3", p, r, s, z)
            else:
                st['pole'] += 1
                for lam, (c1, c2) in PAIRS.items():
                    if vl[lam] == n:
                        st['eq_n'] += 1
                        if max(f[c1].valuation(p), f[c2].valuation(p)) >= 2:
                            st['multizero_pole'] += 1
                    # C7
                    pred = (n - vl[lam]) % 2 if vl[lam] < n else 0
                    assert par[c1] == pred and par[c2] == pred, ("C7", p, r, s, z, lam)
                if vr > 0 or vs > 0:
                    st['pr_ps_pole'] += 1
                    assert n % 2 == 0 and S == 0 and U == 0, ("C5", p, r, s, z)       # C5
                elif n % 2:
                    assert U == 1 and ((vpl > 0) != (vmi > 0)) and S == (1 if vmi > 0 else 0), \
                        ("C6odd", p, r, s, z)
                else:
                    assert S == 0 and U == 0, ("C6even", p, r, s, z)
            # C4 (главное утверждение Шага 2)
            if vT:
                st['S1'] += 1
                assert (r - s) % p == 0, ("C4T", p, r, s, z)
                if (r + s) % p:
                    st['N1_viol'] += 1
            if vL:
                st['SU1'] += 1
                assert (r + s) % p == 0, ("C4L", p, r, s, z)
            # C8
            if p % 4 == 3:
                assert all(x == 0 for x in par.values()), ("C8", p, r, s, z)
            if p == 3:
                assert ez >= 1, ("C8p3", p, r, s, z)
        if idx % 100 == 0:
            print("  p=%d pair %d/%d %s %.0fs" % (p, idx, len(prs), st, time.time() - t0), flush=True)
    print("p=%d DONE pairs=%d %s %.0fs" % (p, len(prs), st, time.time() - t0), flush=True)
    return st


if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(',')]
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
    tot = {}
    for p in ps:
        st = run(p, tl)
        for k, v in st.items():
            tot[k] = tot.get(k, 0) + v
    print("TOTAL", ps, tot, flush=True)
    print("ALL C1-C8 ASSERTIONS PASSED; negative controls: N1_viol=%d (must be >0 for p=1 mod 4), "
          "N2_viol=%d (>0 => the two edge rows are genuinely needed)" % (tot['N1_viol'], tot['N2_viol']), flush=True)
