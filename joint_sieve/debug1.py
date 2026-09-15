from jsieve import *
r, s = 247, 825
z = QQ(168) / QQ(425)**2
cells6 = [s, -s, -r, -(s + r), s - r, -(s - r)]
roots = {lam: (1+lam*z).sqrt() for lam in cells6}
f = Factor([-1072,-825,-578], verbose=True)
P = f.point(z, prod(roots[t] for t in f.T)); print('P', P, 'z', f.z_of(P))
cls = f.decompose(P); print('класс', cls, 'gens', f.gens, 'tgens', f.tgens)
JS = JointSieve(r, s, cells6, [f])
for ell in [7, 23, 29, 31, 37]:
    zb, roots_l, infr = JS.zdata(ell)
    loc = f.local(ell, zb, roots_l, infr)
    El = f.E.change_ring(GF(ell)); G = El.abelian_group()
    Pl = El(P); lp = G.discrete_log(Pl)
    comb = [sum(c*g[k] for c, g in zip(cls, loc['glog'])) % loc['inv'][k] for k in range(2)]
    zl = GF(ell)(z)
    print(ell, 'inv', loc['inv'], 'log P̄', lp, 'комбинация', comb, 'z̄', zl, 'z̄∈Z?', int(zl) in zb, 'fib', loc['fib'].get(int(zl)))
