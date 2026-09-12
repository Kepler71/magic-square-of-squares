# Проверка обратной конструкции Codex (§4) на случайных точках над F_l.
# Для классов (5,13,65) и (6,26,39): случайное p (q=1), корни r,s,w,v,alpha,beta в F_l,
# конструкция t, D, Y1..Y4 -> проверка Y1^2=a+b, Y2^2=a-b, Y3^2=c+b, Y4^2=c-b и (Y3-Y4)=X(Y1-Y2).
import random, functools
print = functools.partial(print, flush=True)
ell = next_prime(10^12); F = GF(ell)
random.seed(int(7))
def trial(d1, d2, d3):
    delta = sqrt(QQ(d1*d2)/d3); assert delta in QQ
    q = F(1)
    while True:
        p = F.random_element()
        vals = {'r': (p^2 - q^2)/d1, 's': (p^2 + q^2)/d2, 'w': (109*p^2 - 229*q^2)/d3, 'v': (229*p^2 - 109*q^2)/d3}
        if not all(x.is_square() and x != 0 for x in vals.values()): continue
        r, s, w, v = [vals[k].sqrt() * random.choice([1, -1]) for k in 'rswv']
        assert 109*q^2 - 229*p^2 == -d3*v^2                       # 4-я квадрика в конвенции Codex
        h = d3*p*q*w*v / (169*(q^4 - p^4))
        if not ((1 + 4*h).is_square() and (1 - 4*h).is_square()): continue
        al = (1 + 4*h).sqrt() * random.choice([1, -1]); be = (1 - 4*h).sqrt() * random.choice([1, -1])
        A = (al + be)/2; B = (al - be)/2
        if A == 1: continue
        t = B/(1 - A)
        D = 2*(t^2 + 1)*q*v/(F(delta)*r*s)
        b = 338*t*(t^2 - 1); X = p/q
        a = 109*(t^2 + 1)^2; c = 229*(t^2 + 1)^2
        Y1 = (2*b/D + D)/2; Y2 = (2*b/D - D)/2
        Y3 = (2*b/(X*D) + X*D)/2; Y4 = (2*b/(X*D) - X*D)/2
        return [Y1^2 == a + b, Y2^2 == a - b, Y3^2 == c + b, Y4^2 == c - b, Y3 - Y4 == X*(Y1 - Y2)]
for cls in [(5, 13, 65), (6, 26, 39)]:
    results = [trial(*cls) for _ in range(200)]
    ok = sum(all(x) for x in results)
    # если не все прошли — какие уравнения падают
    fails = [i for x in results for i, y in enumerate(x) if not y]
    print(f"class {cls}: construction valid in {ok}/200 random trials over F_l;  failing eq. indices: {sorted(set(fails))}")
