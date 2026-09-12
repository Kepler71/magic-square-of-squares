# Сопоставление: Галуа-орбита (char 0) <-> высоты её сечений (mod p).
import functools, sys
print = functools.partial(print, flush=True)
p = int(sys.argv[1]); F = GF(p)
res = load('sections_charpoly.sobj'); fl = res['l1']
Hd = load(f'heights_{p}.sobj'); V = Hd['V']; hs = Hd['hs']
B = V[0][list(V[0].keys())[0]].parent()
Bq = PolynomialRing(QQ, 'a,b,c,d,e,f,g')
a,b,c,d,e,f,g = Bq.gens()
l1 = a + 3*b - 2*c + 5*d + 7*e - 11*f + 13*g
Fx.<T> = F[]
facs = [(h, m) for h, m in fl]
red = [Fx(h.change_ring(F)) for h, _ in facs]
for i in range(len(red)):
    for j in range(i+1, len(red)):
        assert red[i].gcd(red[j]) == 1, "factors collide mod p"
table = {}
for v, hv in zip(V, hs):
    val = l1.change_ring(F)(*[v[x] for x in v[x].parent().gens()]) if False else None
    names = [str(x) for x in B.gens()]
    vals = [v[B(n)] for n in ['a','b','c','d','e','f','g']]
    val = sum(cf*x for cf, x in zip([1,3,-2,5,7,-11,13], vals))
    idx = [i for i, r in enumerate(red) if r(val) == 0]
    assert len(idx) == 1
    table.setdefault(idx[0], []).append(hv)
for i, (h, m) in enumerate(facs):
    fld = "Q" if h.degree() == 1 else str(sorted(Fq.discriminant() for Fq, _, _ in NumberField(h, 'z').subfields() if Fq.degree() == 2))
    print(f"orbit size {h.degree()}, mult {m:2d}, quadratic subfields {fld:28s} heights {sorted(table.get(i, []))}")
