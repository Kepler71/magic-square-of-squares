# Fable, 25.09.2026. Закрытие наклона: полная группа (mwrank proof=True + насыщение max_prime=-1), интервальная матрица Грама,
# лемма ранга 2, случай A (решето по простым + точная проверка), случай B (торсионные исключения). Выход: closure_<p>_<q>.json.
import sys, json, time
sys.path.insert(0, '/home/kep/magicKube/rigorous_fable')
from fmodels import *

def close(fn_cert, sieve_primes=10, verbose=True):
    t0 = time.time(); Dt = json.load(open(fn_cert)); p, q = map(ZZ, Dt['slope'].split('/'))
    Ts = [[QQ(c) for c in T] for T in Dt['Ts']]; hub = Dt['hub']; others = Dt['others']
    Et = EllipticCurve(QQ, list(map(QQ, Dt['E']))); facs = [Fac(T) for T in Ts]; ims = [Image(f, Et) for f in facs]
    S = cells(p, q)
    gens = Et.gens(proof=True, use_database=False); certain = Et.gens_certain()
    print('gens', gens, 'certain', certain, 'secs', time.time() - t0, flush=True); assert certain and len(gens) == 2
    sat, sat_index, _ = Et.saturation(gens, max_prime=-1); assert sat_index == 1
    print('saturation index', sat_index, 'secs', time.time() - t0, flush=True)
    Tor = Et.torsion_points(); ntor = len(Tor)
    I = RealIntervalField(200); Bv = ZZ(Dt['B']); err = I(Bv).log()/3
    def hint(P, n=6):
        Q = (2**n)*P
        if Q.is_zero(): return I(0)
        return I(max(abs(Q[0].numerator()), Q[0].denominator())).log()/4**n + I(-1, 1)*err/4**n
    h1, h2 = hint(gens[0]), hint(gens[1]); ip = (hint(gens[0] + gens[1]) - h1 - h2)/2
    det = h1*h2 - ip*ip; tr = h1 + h2; assert det.lower() > 0
    lam = QQ(floor((det/tr).lower()*64))/64; assert lam > 0 and I(lam).upper() < (det/tr).lower()
    print('Gram', h1, h2, ip, 'lambda', lam, flush=True)
    # константы сравнения высот: |hhat(Q_i) - h(z)| <= C_i
    Cs = {j: ZZ((err + I(ims[j].mobius_norm()).log()).upper().ceil()) for j in range(len(Ts))}
    cpair = {c['j']: QQ(c['integer_bound']) for c in Dt['pairs']}
    j2, j0 = others; c2, c0 = cpair[j2], cpair[j0]; C1 = Cs[hub]; C2, C0 = Cs[j2], Cs[j0]
    print('height constants', Cs, 'pairing constants', {j2: c2, j0: c0}, flush=True)
    def kval(H):
        am = max(c2, c0)**2/(H - C1)
        if H - max(C0, C2) - am <= 0: return None
        return (c2 + c0)**2/(H - C1) + (C0 + C2 + am)**2/(4*(H - max(C0, C2) - am))
    H0 = ZZ(max(Cs.values()) + 1)
    while kval(H0) is None or kval(H0) >= lam: H0 += 1
    bound = ZZ(floor((H0 + C1)/lam)).isqrt() + 1
    print('H0', H0, 'K(H0)', float(kval(H0)), 'box', bound, 'combinations', ntor*(2*bound + 1)**2, flush=True)
    # ---- случай A: H < H0 ⇒ Q_hub = nG1 + mG2 + T, |n|,|m| <= bound
    def zpoint(j, P):
        return ims[j].z_of_x(None if P.is_zero() else P[0])
    def qualifies(z):
        return z is not None and z != 0 and all(1 + k*z > 0 and (1 + k*z).is_square() for k in S)
    triples = [(n, m, k) for n in range(-bound, bound + 1) for m in range(-bound, bound + 1) for k in range(ntor)]
    count = len(triples); sieve_log = []; M = ims[hub].M; G1, G2 = gens
    for pr in prime_range(3, 200):
        if Et.discriminant() % pr == 0 or M.det() % pr == 0: continue
        if any(c.denominator() % pr == 0 for P in gens for c in (P[0], P[1])): continue   # образующая приводится в O: не рискуем
        Fp = GF(pr); Ep = Et.change_ring(Fp); Mp = M.change_ring(Fp)
        def allowed(P):
            X, Z = (Fp(1), Fp(0)) if P.is_zero() else (P[0], Fp(1))
            num = Mp[1,1]*X - Mp[0,1]*Z; den = -Mp[1,0]*X + Mp[0,0]*Z
            if den == 0: return True
            zp = num/den
            return all((1 + Fp(c)*zp).is_square() for c in S)
        good = set(P for P in Ep if allowed(P))
        rp1, rp2 = Ep([Fp(c) for c in G1]), Ep([Fp(c) for c in G2]); rt = [Ep([Fp(c) for c in T]) for T in Tor]
        a = {n: n*rp1 for n in range(-bound, bound + 1)}; b = {(m, k): m*rp2 + rt[k] for m in range(-bound, bound + 1) for k in range(ntor)}
        triples = [(n, m, k) for n, m, k in triples if a[n] + b[m, k] in good]
        sieve_log.append({'p': int(pr), 'allowed_points': len(good), 'order': int(Ep.order()), 'remaining': len(triples)})
        print('sieve p', pr, '|E(Fp)|', Ep.order(), 'allowed', len(good), 'remaining', len(triples), flush=True)
        if len(sieve_log) >= sieve_primes: break
    solA = []; exact = []
    for n, m, k in triples:
        P = n*G1 + m*G2 + Tor[k]; z = zpoint(hub, P)
        if qualifies(z): solA.append(str(z))
        exact.append({'indices': [int(n), int(m), int(k)], 'z': str(z)})
    print('case A: combinations', count, 'exact checks', len(triples), 'solutions', solA, 'secs', time.time() - t0, flush=True)
    # ---- случай B: Q_{j2} - eps Q_{j0} = T ∈ E(Q)_tors ⇒ x(Q_{j0}) = x(Q_{j2} - T)
    Rz = PolynomialRing(QQ, 'z'); Fz = Rz.fraction_field(); z = Fz.gen()
    x2, yc2 = ims[j2].xy(z); x0, _ = ims[j0].xy(z); rad2 = prod(1 + c*z for c in ims[j2].fac.T)
    a2, a4, a6 = Et.a2(), Et.a4(), Et.a6(); f = lambda x: x**3 + a2*x**2 + a4*x + a6
    caseB = []; polys = []
    for T in Tor:
        if T.is_zero(): eqn = (x0 - x2).numerator()
        else:
            tx, ty = T[0], T[1]
            # x(P - T), P=(x,y): lambda = (y + ty)/(x - tx); x(P-T) = lambda^2 - a2 - x - tx = AT + BT*y
            AT = (f(x2) + ty**2)/(x2 - tx)**2 - a2 - x2 - tx; BT = 2*ty/(x2 - tx)**2
            eqn = ((x0 - AT)**2 - BT**2*yc2**2*rad2).numerator()
        assert eqn != 0; polys.append(str(eqn))
        for zz, _ in eqn.roots(QQ): caseB.append(zz)
        # особые точки формул: Q_{j2} = O, Q_{j2} = ±T
        caseB.append(ims[j2].z_of_x(None))
        if not T.is_zero(): caseB.append(ims[j2].z_of_x(T[0]))
    caseB.append(ims[j0].z_of_x(None))
    solB = sorted(set(str(zz) for zz in caseB if qualifies(zz)))
    print('case B candidates', sorted(set(map(str, caseB))), 'solutions', solB, flush=True)
    out = {'slope': Dt['slope'], 'certificate': fn_cert, 'E': Dt['E'], 'gens': [[str(c) for c in P] for P in gens], 'gens_certain': bool(certain),
           'full_saturation_index': int(sat_index), 'torsion': [[str(c) for c in P] for P in Tor],
           'gram_intervals': [[str(QQ(floor(v.lower()*10**8))/10**8), str(QQ(ceil(v.upper()*10**8))/10**8)] for v in (h1, h2, ip)],
           'lambda_lower': str(lam), 'B': str(Bv), 'height_constants': {str(j): int(v) for j, v in Cs.items()}, 'hub': hub, 'others': others,
           'pairing_constants': {str(j2): int(c2), str(j0): int(c0)}, 'H0': int(H0), 'K_H0': str(kval(H0)), 'box_bound': int(bound),
           'case_A_combinations': count, 'sieve': sieve_log, 'exact_candidates': exact, 'case_A_solutions': solA,
           'case_B_polynomials': polys, 'case_B_candidates': sorted(set(map(str, caseB))), 'case_B_solutions': solB, 'seconds': time.time() - t0}
    fn = fn_cert.replace('certificate_', 'closure_'); json.dump(out, open(fn, 'w'), indent=1); print('written', fn, flush=True)
    return out

if __name__ == '__main__':
    for fn in sys.argv[1:]: close(fn)
