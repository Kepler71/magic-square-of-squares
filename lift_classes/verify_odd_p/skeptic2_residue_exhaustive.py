#!/usr/bin/env python3
"""Скептик-2: ПОЛНЫЙ перебор по вычетам (не выборка) для Шага 2 теоремы Codex
об остаточных классах (нечётное p).

Идея. Все локальные условия при нечётном p зависят только от (A,B) = (r z, s z) p-адически.
  * Полюс n = -v_p(z) >= 1:  z = u p^{-n}, u единица; A = r u, B = s u  (не обе делятся на p,
    т.к. gcd(r,s)=1).  p^n f_ij = p^n + i A + j B.  p | (r-s) <=> A = B (mod p), p | (r+s) <=> A = -B.
    Любая пара вычетов (A,B) mod p^M, не обе = 0 mod p, реализуется целыми взаимно простыми r != +-s.
  * z p-целое: A = r z, B = s z  — любая пара (A,B) mod p^M;  f_ij = 1 + i A + j B.
Перебираются ВСЕ (A,B) mod p^M (numpy), n = 1..M-1.  Конфигурации, где какая-то клетка = 0 mod p^M
(оценка не определена на этой точности), пропускаются и считаются отдельно.

Проверки на множестве, где выполнены 8 соотношений ЧЁТНОСТИ (это ровно то, что использует Шаг 2):
  P1 противоположные клетки равной чётности, вектор = матрица (S,U,S+U);
  P2 z целое => все чётности 0;
  P3 чётность v(T) (T считается как сумма оценок трёх клеток) = S, чётность v(L) = S+U;
  P4 v(T) нечётна => A = B (mod p) [p|(r-s)];  v(L) нечётна => A = -B (mod p) [p|(r+s)];
  P5 полюс, p|r или p|s => n чётно и S=U=0;
  P6 полюс, r,s единицы, n нечётно => U=1 и ровно одна из r+-s делится на p; n чётно => S=U=0.
На множестве, где 8 произведений — КВАДРАТЫ в Q_p (полные классы), дополнительно:
  Q1 противоположные клетки одного полного класса;
  Q2 p = 3 mod 4 => все чётности 0 (локальная форма Шага 5).
Отрицательные контроли (обязаны нарушаться, иначе проверка пуста):
  N1 «v(T) нечётна => A = -B» (перепутанная связь) — должно быть нарушено при p = 1 mod 4;
  N2 только 4 центральные линии + 2 строки (без столбцов) — P4 должно нарушаться.
Все циклы конечны (по p, n, k < M).  Аргументы: p:M через запятую.
"""
import sys, time
import numpy as np

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IX = {c: k for k, c in enumerate(CELLS)}
ROWS = [[IX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)]
COLS = [[IX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)]
DIAG = [[IX[(-1, -1)], IX[(0, 0)], IX[(1, 1)]], [IX[(-1, 1)], IX[(0, 0)], IX[(1, -1)]]]
LINES = ROWS + COLS + DIAG
LINES_N2 = ROWS + DIAG + [COLS[1]]          # без двух крайних столбцов
T_CELLS = [IX[(-1, 0)], IX[(1, 1)], IX[(0, -1)]]
L_CELLS = [IX[(-1, 0)], IX[(1, -1)], IX[(0, 1)]]


def val_and_chi(X, p, M, chi):
    """X: int64 массив вычетов mod p^M. Возвращает (v, chi(unit part)), v=M означает 'не определено'."""
    v = np.zeros(X.shape, dtype=np.int8)
    cur = X.copy()
    alive = np.ones(X.shape, dtype=bool)
    for _k in range(M):                      # явная граница
        d = alive & (cur % p == 0)
        v[d] += 1
        cur[d] //= p
        alive = d
    c = chi[(cur % p).astype(np.int64)]
    return v, c


def run(p, M):
    t0 = time.time()
    q = p ** M
    chi = np.array([0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)], dtype=np.int8)
    A = np.arange(q, dtype=np.int64)[:, None] * np.ones((1, q), dtype=np.int64)
    B = np.arange(q, dtype=np.int64)[None, :] * np.ones((q, 1), dtype=np.int64)
    Ap, Bp = A % p, B % p
    rs_eq = (Ap == Bp)            # p | (r-s)
    rs_op = ((Ap + Bp) % p == 0)  # p | (r+s)
    stats = {}
    # режимы: n=0 означает z p-целое (f = 1 + iA + jB); n>=1 — полюс
    for n in range(0, M):
        base = 1 if n == 0 else p ** n
        V = []; C = []
        for (i, j) in CELLS:
            X = (base + i * A + j * B) % q
            v, c = val_and_chi(X, p, M, chi)
            V.append(v); C.append(c)
        det = np.ones(A.shape, dtype=bool)
        for v in V:
            det &= (v < M)
        if n >= 1:
            det &= ~((Ap == 0) & (Bp == 0))     # взаимная простота r,s
        par = [((v.astype(np.int16) - n) % 2).astype(np.int8) for v in V]
        linepar = np.ones(A.shape, dtype=bool)
        for Ln in LINES:
            linepar &= ((par[Ln[0]] + par[Ln[1]] + par[Ln[2]]) % 2 == 0)
        linesq = linepar.copy()
        for Ln in LINES:
            linesq &= (C[Ln[0]] * C[Ln[1]] * C[Ln[2]] == 1)
        mp = det & linepar
        ms = det & linesq
        S = par[IX[(1, 1)]]; U = par[IX[(1, 0)]]; SU = par[IX[(1, -1)]]
        # P1
        for (i, j) in CELLS:
            assert not np.any(mp & (par[IX[(i, j)]] != par[IX[(-i, -j)]])), ("P1 opp", p, n, i, j)
        assert not np.any(mp & (SU != (S + U) % 2)), ("P1 matrix", p, n)
        assert not np.any(mp & (par[IX[(0, 1)]] != U)), ("P1 U", p, n)
        anyodd = np.zeros(A.shape, dtype=bool)
        for pr in par:
            anyodd |= (pr == 1)
        # P2
        if n == 0:
            assert not np.any(mp & anyodd), ("P2 integral odd", p)
        vT = (sum(V[k].astype(np.int16) for k in T_CELLS) - 3 * n) % 2
        vL = (sum(V[k].astype(np.int16) for k in L_CELLS) - 3 * n) % 2
        # P3
        assert not np.any(mp & (vT != S)), ("P3 T", p, n)
        assert not np.any(mp & (vL != SU)), ("P3 L", p, n)
        # P4
        assert not np.any(mp & (vT == 1) & ~rs_eq), ("P4 T", p, n)
        assert not np.any(mp & (vL == 1) & ~rs_op), ("P4 L", p, n)
        if n >= 1:
            pr_ps = (Ap == 0) | (Bp == 0)
            # P5
            if n % 2 == 1:
                assert not np.any(mp & pr_ps), ("P5 n odd with p|r or p|s", p, n)
            assert not np.any(mp & pr_ps & anyodd), ("P5 S=U=0", p, n)
            # P6
            units = ~pr_ps
            if n % 2 == 1:
                assert not np.any(mp & units & (U != 1)), ("P6 U", p, n)
                assert not np.any(mp & units & ~(rs_eq ^ rs_op)), ("P6 exactly one", p, n)
            else:
                assert not np.any(mp & units & anyodd), ("P6 even n", p, n)
        # Q1: полный класс противоположных клеток
        for (i, j) in CELLS:
            a, b = IX[(i, j)], IX[(-i, -j)]
            assert not np.any(ms & (C[a] != C[b])), ("Q1", p, n, i, j)
        # Q2
        if p % 4 == 3:
            assert not np.any(ms & anyodd), ("Q2", p, n)
        # отрицательные контроли
        n1 = int(np.sum(ms & (vT == 1) & ~rs_op))
        linep2 = np.ones(A.shape, dtype=bool)
        for Ln in LINES_N2:
            linep2 &= ((par[Ln[0]] + par[Ln[1]] + par[Ln[2]]) % 2 == 0)
        n2 = int(np.sum(det & linep2 & (((vT == 1) & ~rs_eq) | ((vL == 1) & ~rs_op))))
        # покрытие: кратные нули и случай v(lambda)=n
        multi = np.zeros(A.shape, dtype=bool)
        for v in V:
            multi |= (v.astype(np.int16) - n >= 2) if n == 0 else (v.astype(np.int16) >= n + 2)
        st = dict(n=n, total=int(A.size), undetermined=int(np.sum(~det)), par_ok=int(mp.sum()),
                  sq_ok=int(ms.sum()), par_nontriv=int(np.sum(mp & anyodd)),
                  sq_nontriv=int(np.sum(ms & anyodd)),
                  sq_T_odd=int(np.sum(ms & (vT == 1))), sq_L_odd=int(np.sum(ms & (vL == 1))),
                  par_multizero=int(np.sum(mp & multi)), sq_multizero=int(np.sum(ms & multi)),
                  N1_viol=n1, N2_viol=n2)
        if n >= 1:
            st['par_pole_pr_ps'] = int(np.sum(mp & ((Ap == 0) | (Bp == 0))))
        stats[n] = st
        print(f"  p={p} M={M} n={n}: {st}  {time.time()-t0:.1f}s", flush=True)
    return stats


if __name__ == '__main__':
    spec = sys.argv[1] if len(sys.argv) > 1 else "3:5,5:3"
    allst = {}
    for item in spec.split(','):
        p, M = map(int, item.split(':'))
        assert p % 2 == 1 and p ** (2 * M) <= 7 * 10 ** 6, "слишком большой перебор"
        print(f"p={p}, M={M}: pairs={(p**M)**2}", flush=True)
        allst[p] = run(p, M)
        tot_nt = sum(s['sq_nontriv'] for s in allst[p].values())
        n1 = sum(s['N1_viol'] for s in allst[p].values())
        n2 = sum(s['N2_viol'] for s in allst[p].values())
        print(f"p={p} DONE: Q_p-nontrivial={tot_nt}, N1_viol={n1}, N2_viol={n2}", flush=True)
    print("ALL ASSERTIONS PASSED", flush=True)
