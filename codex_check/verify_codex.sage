# Независимая проверка NOTE_FOR_CLAUDE_2026-09-11_FROM_CODEX.md, §1–2 (Claude Code).
import itertools, functools
print = functools.partial(print, flush=True)

# ---- (1) новая точка на Q: Y^2 = (109X^2-229)(X^4-1) и на квадриках класса (-15,2,-30)
X0 = QQ(26399)/54721; Y0 = QQ(2274007548365280)/163855897047361
print("(1) (X0,Y0) on Q:", Y0^2 == (109*X0^2 - 229)*(X0^4 - 1))
p, q, r, sc, w = 26399, 54721, 12376, 42961, 142566
print("    p^2-q^2 = -15r^2:", p^2 - q^2 == -15*r^2,
      "| p^2+q^2 = 2s^2:", p^2 + q^2 == 2*sc^2,
      "| 109p^2-229q^2 = -30w^2:", 109*p^2 - 229*q^2 == -30*w^2, "| gcd(p,q) =", gcd(p, q))

# ---- (2) тождество исключения, выведено заново:
# D = Y1-Y2, u = D^2; из Y1^2 = a+b, Y2^2 = a-b:  u^2 - 4au + 4b^2 = 0
# из Y3^2 = c+b, Y4^2 = c-b, Y3-Y4 = X D:        X^4 u^2 - 4cX^2 u + 4b^2 = 0
R.<t, X, u> = QQ[]
a = 109*(t^2+1)^2; b = 338*t*(t^2-1); c = 229*(t^2+1)^2
e1 = u^2 - 4*a*u + 4*b^2
e2 = X^4*u^2 - 4*c*X^2*u + 4*b^2
mine = (b*(X^4 - 1))^2 - 4*X^2*(c*X^2 - a)*(a*X^2 - c)       # мой результат исключения
codex = (338*t*(t^2-1)*(1 - X^4))^2 - (229 - 109*X^2)*(109 - 229*X^2)*(2*X*(t^2+1)^2)^2
print("(2) Codex identity == mine (as polynomials):", codex == mine)
res = e1.resultant(e2, u)
print("    mine divides Res_u(e1,e2):", mine.divides(res), "  Res factors:", [(f.degree(), m) for f, m in res.factor()])
# инвариантность условия относительно X -> 1/X (конвенция p/q vs q/p не важна)
S.<x> = QQ[]
g4 = (229 - 109*x^2)*(109 - 229*x^2)
print("    x^4 * g4(1/x) == g4(x):", S(x^4 * g4(1/x)) == g4)

# ---- (3) запрет по модулю 7 для дополненной системы класса (-15,2,-30)
def eqs(P, Q_, Rr, Sc, W, V, m):
    return [(P^2 - Q_^2 + 15*Rr^2) % m, (P^2 + Q_^2 - 2*Sc^2) % m,
            (109*P^2 - 229*Q_^2 + 30*W^2) % m, (109*Q_^2 - 229*P^2 - 30*V^2) % m]
sol3 = sol4 = 0; q0 = 0
for P, Q_, Rr, Sc, W, V in itertools.product(range(7), repeat=6):
    if (P, Q_, Rr, Sc, W, V) == (0,)*6: continue
    ev = eqs(P, Q_, Rr, Sc, W, V, 7)
    if ev[:3] == [0, 0, 0]:
        sol3 += 1
        if ev[3] == 0: sol4 += 1
print("(3) nonzero solutions mod 7: first three quadrics:", sol3, "(контроль: должно быть > 0)",
      "| all four:", sol4, "(Codex: 0)")
# то же для класса (15,2,30): p<->q, w<->v
sol4b = 0
for P, Q_, Rr, Sc, W, V in itertools.product(range(7), repeat=6):
    if (P, Q_, Rr, Sc, W, V) == (0,)*6: continue
    if [(P^2 - Q_^2 - 15*Rr^2) % 7, (P^2 + Q_^2 - 2*Sc^2) % 7,
        (109*P^2 - 229*Q_^2 - 30*W^2) % 7, (109*Q_^2 - 229*P^2 + 30*V^2) % 7] == [0]*4:
        sol4b += 1
print("    class (15,2,30), all four quadrics, nonzero solutions mod 7:", sol4b)
# контроль: новая точка удовлетворяет первым трём mod 7 ненулевым образом
print("    new point mod 7 (first three):", eqs(p, q, r, sc, w, 0, 7)[:3], " reduction nonzero:", any(x % 7 for x in (p, q, r, sc, w)))

# ---- (4) согласованность: новая точка НЕ должна удовлетворять 4-й квадрике
val = (109*q^2 - 229*p^2) / QQ(30)
print("(4) (109q^2-229p^2)/30 is a rational square:", val.is_square(), "(ожидается False)")
print("    (229q^2-109p^2)(109q^2-229p^2) is a square:", ((229*q^2 - 109*p^2)*(109*q^2 - 229*p^2)).is_square())
