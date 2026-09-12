# 2-спуск над квадратичным полем k для E с ОДНОЙ рациональной точкой порядка 2 (Claude, 2026-09-12).
# E: y^2 = x(x^2 + a x + b), Delta = a^2 - 4b не квадрат в k. Этальная алгебра L = k x K, K = k(sqrt Delta) = k[t]/(t^2+a t+b).
# mu: E(k)/2E(k) -> L*/L*^2,  P=(x,y) -> (x, x - theta);  для T=(0,0): (b, -theta)  (правило f'(r) при x = r).
# Образ лежит в ядре нормы N_{L/k}: (alpha, beta) -> alpha * N_{K/k}(beta) в k*/k*^2.
# Локально: dim E(k_v)/2E(k_v) = dim E(k_v)[2] + [k_v:Q_2] (v | 2), = dim E(k_v)[2] (конечные v ∤ 2),
#           = dim E(R)[2] - 1 (вещественные). Локальные образы — выборка до точной размерности (иначе ошибка).
# Граница: rank E(k) <= dim Sel^2 - dim E(k)[2] = dim Sel^2 - 1.
import random, functools
print = functools.partial(print, flush=True)

class PartialDescent:
    def __init__(self, k, a, b, extra_primes=()):
        self.k = k; self.a = k(a); self.b = k(b)
        Rx = PolynomialRing(k, 'X'); X = Rx.gen(); self.Rx = Rx
        self.g = X^2 + self.a*X + self.b
        self.f = X*self.g
        assert not (self.a^2 - 4*self.b).is_square(), "Delta — квадрат: применим полный спуск"
        self.K = k.extension(self.g, 'th'); self.th = self.K.gen()
        self.Kabs = self.K.absolute_field('ta'); self.toKabs, self.fromKabs = self.Kabs.structure()[1], self.Kabs.structure()[0]
        disc = self.b * (self.a^2 - 4*self.b)
        ps = set(ZZ(2*disc.norm()).prime_factors()) | set(extra_primes)
        self.S = sorted(sum([k.primes_above(p) for p in sorted(ps)], []), key=lambda P: (P.smallest_integer(), str(P)))
        self.SK = sorted(sum([self.Kabs.primes_above(p) for p in sorted(ps)], []), key=lambda P: (P.smallest_integer(), str(P)))
        self.gk = k.selmer_generators(self.S, 2)
        self.gK = self.Kabs.selmer_generators(self.SK, 2)
        self.nk, self.nK = len(self.gk), len(self.gK)
        self.plk = [LocalSq(k, P) for P in self.S]
        self.emb_k = k.real_embeddings()

    # ---- координаты локальных классов ----
    def kcoords(self, alpha, i):
        return self.plk[i].coords(alpha)
    def Kcoords_at(self, beta, P):
        """ координаты beta в ⊕_{W | P} K_W*/K_W*^2 """
        out = []
        for W in self.Kabs.primes_above(P.smallest_integer()):
            if not self._above(W, P): continue
            out += LocalSq(self.Kabs, W).coords(self.toKabs(beta) if beta.parent() is self.K else self.Kabs(beta))
        return out
    def _above(self, W, P):
        if not hasattr(self, '_abv'): self._abv = {}
        key = (str(W), str(P))
        if key not in self._abv:
            # W лежит над P, если образы образующих P лежат в W
            self._abv[key] = all(self.toKabs(self.K(self.k(g))).valuation(W) > 0 for g in P.gens() if g != 0)
        return self._abv[key]
