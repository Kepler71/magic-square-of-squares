# Claude, 25.09: независимый аудит закрытия 126/451 (Codex). Код Codex не импортируется.
from sage.all import *
import itertools, json, time
r, s = 126, 451
S = [s, -s, r, -r, s - r, r - s, s + r, -(s + r)]      # клетки 1+c z
E = EllipticCurve([0, 0, 0, -21739732, 23939545056])
T = {0: [-577, -451, -325, -126], 1: [-325, -126, 126], 2: [-126, 126, 325]}
M = {0: [[-2730602, 7490], [-577, 1]], 1: [[-3439800, 100333], [-2925, 9]], 2: [[647388, 7490], [-126, 1]]}
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Kz = Rz.fraction_field()
def mob(m, t): return (m[0][0]*t + m[0][1]) / (m[1][0]*t + m[1][1])
out = {}
# (а) отображения: f(M_i(z)) / prod(1+cz) — квадрат в Q(z)  (т.е. y = sqrt(prod)*rational)
fE = lambda X: X**3 - 21739732*X + 23939545056
for i in (0, 1, 2):
    X = mob(M[i], Kz(z))
    ratio = Kz(fE(X)) / Kz(prod(1 + c*z for c in T[i]))
    num, den = ratio.numerator(), ratio.denominator()
    ok = num.is_square() and den.is_square() if num.leading_coefficient() > 0 else False
    # допускаем постоянный множитель-квадрат в Q
    lc = QQ(num.leading_coefficient()) / QQ(den.leading_coefficient())
    ok = (Rz(num / num.leading_coefficient()).is_square() and Rz(den / den.leading_coefficient()).is_square() and lc.is_square())
    out[f'map{i}'] = bool(ok)
print('(а) отображения z -> E корректны:', out, flush=True)
# (б) ранг и образующие — другим путём: PARI ellrank + Sage saturation
lo, hi, *_ = [int(t) if k < 2 else t for k, t in enumerate(pari(E).ellrank())]
G = [E(-2883, 250305), E(-1610, 234024)]
sat, idx, reg = E.saturation(G, max_prime=0)
tors = E.torsion_points()
print(f'(б) PARI ellrank [{lo},{hi}], индекс насыщения {idx}, кручение {len(tors)}, регулятор {float(reg):.4f}', flush=True)
out['rank'] = (lo, hi); out['sat_index'] = int(idx); out['ntors'] = len(tors)
# минимальное собственное число матрицы Грама (численно, контроль к интервальной оценке Codex 1.5)
Gm = matrix(RR, 2, 2, [[G[a].height() if a == b else (( (G[a]+G[b]).height() - (G[a]-G[b]).height()) / 4) for b in range(2)] for a in range(2)])
ev = min(Gm.eigenvalues())
print(f'   матрица Грама (нормировка Sage): {Gm.list()}, мин. собств. {float(ev):.4f}', flush=True)
# (в) K(H) из леммы ранга 2 — точным пересчётом
def Kf(H):
    H = QQ(H); A = QQ(4900) / (H - 45)
    return QQ(9604) / (H - 45) + (87 + A)**2 / (4 * (H - 44 - A))
out['K7728'] = str(Kf(7728)); out['K7728_lt_1.5'] = bool(Kf(7728) < QQ(3)/2)
print(f'(в) K(7728) = {Kf(7728)} < 3/2: {Kf(7728) < QQ(3)/2}', flush=True)
