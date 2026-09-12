import functools
print = functools.partial(print, flush=True)
k2.<r> = QuadraticField(15)
E = EllipticCurve(k2, [0, 0, 0, -100326681, 98571779592*r])
G = E(-676*r - 8619, 30420*r + 182520)
tam = [(P.norm(), E.local_data(P).tamagawa_number(), E.local_data(P).kodaira_symbol()) for P in E.conductor().prime_factors()]
print("bad primes (norm, Tamagawa c_v, Kodaira):", tam)
c = lcm([t_[1] for t_ in tam])
print("exponent bound c (lcm of c_v):", c, "  torsion:", E.torsion_subgroup().invariants())
# точка c*G лежит в E^gr (компонентные группы — порядка c_v, их показатель делит c_v)
cG = c*G
hG, hcG = G.height(), cG.height()
print("h(G) =", hG, "  h(cG) =", hcG, "  ratio =", hcG/hG, "(= c^2 =", c^2, ")")
HF = E.height_function()
lam = HF.min_gr(0.1, 5)
print("lambda_gr (Sage min_gr):", lam)
print("normalization check: h(cG) >= lambda_gr ?", hcG >= lam, "; h(cG)*2 >= lambda_gr ?", 2*hcG >= lam)
# индекс насыщения <= c * sqrt(h(G)/lambda) — берём худшую из двух нормировок (lambda/2) для надёжности
for nm, lam_ in [("same normalization", lam), ("pessimistic lambda/2", lam/2)]:
    B = c*sqrt(hG/lam_)
    print(f"  index bound ({nm}): c*sqrt(h/lambda) = {B}")
B = floor(c*sqrt(hG/(lam/2)))
print("primes to check up to:", B)
bad = [pp for pp in primes(2, B + 1) if E.saturation([G], one_prime=pp)[1] != 1]
print("p-saturation fails at:", bad if bad else "none")
