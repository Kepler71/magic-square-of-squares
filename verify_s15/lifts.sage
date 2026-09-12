# П.3: над X in {±3/2, ±5, ±2/3, ±1/5} все рациональные точки C_{1/5} вырождены (клетки повторяются).
# Независимо от критерия Codex: из необходимого тождества (1)
#   [338 t(t^2-1)(1-X^4) / (2X(t^2+1)^2)]^2 = g4(X) = (229-109X^2)(109-229X^2)
# при фиксированном X получаем t(t^2-1)/(t^2+1)^2 = ±h0, h0 = X sqrt(g4)/(169(1-X^4)); ищем ВСЕ рациональные t,
# затем прямо проверяем все Y_i и равенство X = (Y3-Y4)/(Y1-Y2) при каком-то выборе знаков.
# Клетки (v5-detailed §20.3–20.4): 1, B(s)^2, H(s)^2, B(t)^2, H(t)^2, (H(s)^2+H(t)^2)/2, (H(s)^2+B(t)^2)/2,
# (B(s)^2+H(t)^2)/2, (B(s)^2+B(t)^2)/2, s = 1/5, B(u) = (1+2u-u^2)/(1+u^2), H(u) = (1-2u-u^2)/(1+u^2).
import functools, itertools
print = functools.partial(print, flush=True)
S.<T> = QQ[]
a = 109*(T^2+1)^2; b = 338*T*(T^2-1); c = 229*(T^2+1)^2
# контроль формулы клеток <-> Y_i (v5 §21.1): Y1^2 = 109t^4+338t^3+218t^2-338t+109 = a+b и т.д.
assert a + b == 109*T^4 + 338*T^3 + 218*T^2 - 338*T + 109 and c + b == 229*T^4 + 338*T^3 + 458*T^2 - 338*T + 229
Bf = lambda u: (1 + 2*u - u^2)/(1 + u^2); Hf = lambda u: (1 - 2*u - u^2)/(1 + u^2)
s0 = QQ(1)/5
Fr = S.fraction_field(); tt = Fr(T)
corners = [(Hf(s0)^2 + Hf(tt)^2)/2, (Hf(s0)^2 + Bf(tt)^2)/2, (Bf(s0)^2 + Hf(tt)^2)/2, (Bf(s0)^2 + Bf(tt)^2)/2]
# углы * 169 (1+t^2)^2 должны быть {a±b, c±b} (как множество)
cset = sorted([S((169*(1 + T^2)^2 * cr).numerator()) for cr in corners], key=str)
print("corners*169(1+t^2)^2 == {a±b, c±b}:", sorted([a + b, a - b, c + b, c - b], key=str) == cset)

def cells(t0):
    return [QQ(1), Bf(s0)^2, Hf(s0)^2, Bf(t0)^2, Hf(t0)^2,
            (Hf(s0)^2 + Hf(t0)^2)/2, (Hf(s0)^2 + Bf(t0)^2)/2, (Bf(s0)^2 + Hf(t0)^2)/2, (Bf(s0)^2 + Bf(t0)^2)/2]

Xs = [s_*x for s_ in (1, -1) for x in (QQ(3)/2, QQ(5), QQ(2)/3, QQ(1)/5)]
g4 = lambda x: (229 - 109*x^2)*(109 - 229*x^2)
allt = set()
for X0 in Xs:
    G = g4(X0)
    assert G.is_square(), f"g4({X0}) not a square — X0 would not lift at all"
    h0 = X0*G.sqrt()/(169*(1 - X0^4))
    ts = set()
    for sg in (1, -1):
        quart = T*(T^2 - 1) - sg*h0*(T^2 + 1)^2
        ts |= set(quart.roots(QQ, multiplicities=False))
    lifts = []
    for t0 in sorted(ts):
        vals = [a(t0) + b(t0), a(t0) - b(t0), c(t0) + b(t0), c(t0) - b(t0)]
        if not all(v.is_square() for v in vals): continue
        Y = [v.sqrt() for v in vals]
        ok = False
        for sg in itertools.product((1, -1), repeat=4):
            y = [e*v for e, v in zip(sg, Y)]
            if y[0] != y[1] and (y[2] - y[3])/(y[0] - y[1]) == X0:
                ok = True; break
        if ok:
            cl = cells(t0)
            assert all(x_.is_square() for x_ in cl)
            lifts.append((t0, len(set(cl)), sorted(set(cl))))
            allt.add(t0)
    print(f"X = {str(X0):6s}: rational t from (1): {sorted(ts)}; lifts: {[(l[0], l[1]) for l in lifts]}")
    for l in lifts:
        assert l[1] < 9, f"NON-DEGENERATE LIFT at X={X0}, t={l[0]}"
print("\nall lifting t:", sorted(allt))
print("matches v5 §21.2 list {±1/5, ±5, ±2/3, ±3/2}:", sorted(allt) == sorted([s_*x for s_ in (1, -1) for x in (QQ(1)/5, 5, QQ(2)/3, QQ(3)/2)]))
print("distinct cell values over all lifts:", sorted(set(v for t0 in allt for v in cells(t0))))
