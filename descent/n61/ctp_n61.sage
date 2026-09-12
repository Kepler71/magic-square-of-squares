# CTP на Sel^2(E1/k), k = Q(sqrt165), сечение (71,49,61). rank E1(k) <= dim Sel^2 - 2 - rank(CTP).
import sys, time
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else int(1))
k.<r> = QuadraticField(165)
rts = [2939651, 163724*r, -2939651]
Cv = Curve3(k, rts)
Sel = Cv.selmer_full()
print("dim Sel2 =", Sel.dimension())
els = []
for w in Sel.basis():
    a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n)); b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
    els.append((k(a), k(b)))
# известная точка G (модель Codex, x = 244^2 X): образ в модели Cv
q = 175015061936
Eo = EllipticCurve(k, [0, -9747472064*r, 0, -q^2, 9747472064*r*q^2])
G = Eo(-173017629136, -784866450593280*r - 4316765478263040)
xG = G[0] / 244^2
e1, e2, e3 = Cv.e
assert Cv.f(xG) == (G[1]/244^3)^2
kG = Cv.kappa(xG)
tors = [Cv.kappa(e1), Cv.kappa(e2)]
ct = CTP(Cv); ct.randomize_conic = True; ct.cert = []
t0 = time.time()
n = len(els); M = matrix(GF(2), n, n); CERT_ROWS = []
for i in range(n):
    ct.cert = []
    vals, nP = ct.pair(els[i], els + tors + [kG], reps=2)
    CERT_ROWS.append({'row': i, 'eps': [str(els[i][0]), str(els[i][1])], 'values': vals,
                      'witnesses': ct.witness_head, 'places': ct.cert})
    for j in range(n): M[i, j] = vals[j]
    print(f"row {i}: {vals[:n]} | torsion,G: {vals[n:]} | places {nP} | {time.time()-t0:.0f}s")
print(M)
print("symmetric:", M == M.transpose(), " zero diagonal:", all(M[i, i] == 0 for i in range(n)))
print("rank CTP =", M.rank(), " => rank E1(k) <=", Sel.dimension() - 2 - M.rank())
# независимость образов T1, T2, G: по локальным координатам во всех местах S (достаточное условие)
rows = []
for pr in [tors[0], tors[1], kG]:
    v = []
    for pl in Cv.places: v += pl.coords(pr[0]) + pl.coords(pr[1])
    rows.append(v)
print("local rank of images of T1, T2, G:", matrix(GF(2), rows).rank(), "(3 => independent => kernel of CTP = <T1, T2, G>)")

# свидетельства локальных образов полного 2-спуска (представители и достигнутые размерности)
SEL_WITNESS = []
for pl in Cv.places:
    reps = Cv._reps.get(id(pl), [])
    W = matrix(GF(2), [pl.coords(a_) + pl.coords(b_) for a_, b_ in reps])
    SEL_WITNESS.append({'place': str(getattr(pl, 'P', 'real')), 'target_dim': int(Cv.target_full(pl)),
                        'achieved_dim': int(W.rank()), 'image_reps': [[str(a_), str(b_)] for a_, b_ in reps],
                        'local_dim_kv': int(pl.dim)})
import json, hashlib
src = {f: hashlib.sha256(open(f, 'rb').read()).hexdigest() for f in ['/home/kep/magicKube/descent/ek_descent.sage', '/home/kep/magicKube/descent/ctp.sage', '/home/kep/magicKube/descent/n61/ctp_n61.sage']}
cert = {'model_roots': [str(t) for t in rts], 'k': 'Q(sqrt165)', 'S': [str(P) for P in Cv.S],
        'kS2_basis': [str(g) for g in Cv.gens], 'sel2_basis': [[str(a), str(b)] for a, b in els],
        'torsion_images': [[str(a), str(b)] for a, b in tors], 'G_image': [str(kG[0]), str(kG[1])],
        'matrix': [[int(M[i, j]) for j in range(n)] for i in range(n)], 'rank_CTP': int(M.rank()),
        'rank_bound': int(Sel.dimension() - 2 - M.rank()), 'rows': CERT_ROWS, 'sel2_local_witness': SEL_WITNESS, 'sources_sha256': src,
        'seed': int(sys.argv[1]) if len(sys.argv) > 1 else 1}
import os
os.makedirs('/home/kep/magicKube/descent/n61/snapshot', exist_ok=True)
for f in src:
    open('/home/kep/magicKube/descent/n61/snapshot/' + f.split('/')[-1], 'w').write(open(f).read())
json.dump(cert, open(f'/home/kep/magicKube/descent/n61/ctp_n61_certificate_full_seed{cert["seed"]}.json', 'w'), indent=1, default=str)
print("certificate written")
