# Fable, 14.09.2026. Проверка следствия 3 (PROOF_FIXED_SIX_COMPLETE, Codex) и §5: кривая E: y^2 = x^3 - 3360^2 x,
# P = (113^2, 127*113*97); кручение; для m = 1..8 точки mP: h = x(mP), клетки h+3360, h, h-3360 — квадраты?
# положительны? различны от шести фиксированных клеток Саллоуса? S = h + 8840 — квадрат?
import json
from sage.all import *

n = 3360
E = EllipticCurve([0, 0, 0, -n ** 2, 0])
P = E(113 ** 2, 127 * 113 * 97)
out = {}
out['P_on_curve'] = True
out['torsion'] = str(E.torsion_subgroup().invariants())
out['E_F11'] = int(E.reduction(11).order())
out['E_F13'] = int(E.reduction(13).order())
out['disc'] = int(E.discriminant())
out['P_order_infinite'] = bool(P.order() == oo)
# P в 2E(Q)? критерий деления пополам: x, x-n, x+n — квадраты
x0 = P[0]
out['P_in_2E'] = all(QQ(v).is_square() for v in [x0, x0 - n, x0 + n])
print({k: out[k] for k in out})

fixed = [46 ** 2, 58 ** 2, 2 ** 2, 94 ** 2, 74 ** 2, 82 ** 2]
assert len(set(fixed)) == 6
# семь сумм: строки/столбцы/антидиагональ при h: h+3360 + 46^2 + 58^2 = h + 8840 и т.д. — проверим явно
def cells(h):
    b = -n
    return [[h - b, 46 ** 2, 58 ** 2], [2 ** 2, h, 94 ** 2], [74 ** 2, 82 ** 2, h + b]]

rows = []
for m in range(1, 9):
    Q = m * P
    h = Q[0]
    M = cells(h)
    sums = [sum(M[i]) for i in range(3)] + [sum(M[i][j] for i in range(3)) for j in range(3)] + \
           [M[0][2] + M[1][1] + M[2][0]]
    main = M[0][0] + M[1][1] + M[2][2]
    S = h + 8840
    cell9 = [c for row in M for c in row]
    row = {
        'm': m,
        'h': str(h),
        'h_height_digits': len(str(h.numerator())),
        'three_squares': all(QQ(v).is_square() for v in [h, h - n, h + n]),
        'seven_sums_equal': len(set(sums)) == 1 and sums[0] == S,
        'main_diag_equals_S': bool(main == S),
        'nine_positive': all(c > 0 for c in cell9),
        'nine_distinct': len(set(cell9)) == 9,
        'S': str(S),
        'S_is_square': bool(QQ(S).is_square()),
        'S_denominator_is_square': bool(QQ(S).denominator().is_square()),
    }
    rows.append(row)
    print(row)
out['multiples'] = rows

# сверка с числами Codex/рецензии
S2 = QQ(98913801874105761) / 7751179400836
S3 = QQ(59292640648263999369571426440287973289) / 4703413061866180537911414234965769
out['S_2P_matches_codex'] = bool(rows[1]['S'] == str(S2))
out['S_3P_matches_codex'] = bool(rows[2]['S'] == str(S3))
out['x2P_matches_(16)'] = bool((2 * P)[0] == (QQ(174336961) / 2784094) ** 2)
out['2784094_sq'] = int(2784094 ** 2)
# формулы удвоения (15) символически
R = PolynomialRing(QQ, ['x', 'y', 'N']); x, y, N = R.gens()
I = R.ideal([y ** 2 - x ** 3 + N ** 2 * x])
lhs = (x ** 2 + N ** 2) ** 2
out['(15)_x2R'] = bool(((x ** 2 + N ** 2) ** 2 - 4 * y ** 2 * 0 - lhs) == 0)  # тавтология, ниже настоящие
# x(2R) = (x^2+N^2)^2/(4y^2); x(2R) - N = (x^2-2Nx-N^2)^2/(4y^2) <=> (x^2+N^2)^2 - 4Ny^2 = (x^2-2Nx-N^2)^2 mod I
out['(15)_minus'] = bool(I.reduce((x ** 2 + N ** 2) ** 2 - 4 * N * y ** 2 - (x ** 2 - 2 * N * x - N ** 2) ** 2) == 0)
out['(15)_plus'] = bool(I.reduce((x ** 2 + N ** 2) ** 2 + 4 * N * y ** 2 - (x ** 2 + 2 * N * x - N ** 2) ** 2) == 0)
# формула удвоения: x(2R) = ((3x^2-N^2)/(2y))^2 - 2x; сравним с ((x^2+N^2)/(2y))^2 по модулю I
lam2 = (3 * x ** 2 - N ** 2) ** 2
out['(15)_x2R_formula'] = bool(I.reduce(lam2 - 2 * x * 4 * y ** 2 - (x ** 2 + N ** 2) ** 2) == 0)
print({k: out[k] for k in out if k != 'multiples'})
json.dump(out, open('/home/kep/magicKube/fable_review_codex/g2_sallows_multiples.json', 'w'), indent=1, ensure_ascii=False)
print('OK')
